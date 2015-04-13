#! /bin/bash -x

#Add structural variations to reference
./tools/SVSim/SVSim ./data/SVInput.txt ./data/DH10B-K12.fasta
#First you need to do is simulate reads from reference with SV
#For this purpose we used art - http://www.niehs.nih.gov/research/resources/software/biostatistics/art/  
#Please export art_illumina to paths
art_illumina -i ./data/DH10B-K12.sim.fasta -p -f 35 -l 100 -m 300 -s 50 -qs 50 -qs2 50 -o ./data/R
#Then for read correction start Quake
echo -e "./data/R1.fq ./data/R2.fq" >./data/readnames.txt
cat ./data/R1.fq ./data/R2.fq | count-qmers -k 17 -q 64 >counts.txt
#as we know coverage lets chose value like 6, because quake often fails to di it by itself
correct -f ./data/readnames.txt -k 17 -c 6 -m counts.txt -p 4
#Reads are stored in ./data/R1.cor.fq and ./data/R2.cor.fq 
rm kmers.txt error_model* r.log counts.txt

#Now we start assembly with SGA
#----------------ASSEMBLY STEP----------------------------

#
# Example assembly of 100bp E.coli data set. 
#
IN1=./data/R1.cor.fq
IN2=./data/R2.cor.fq

# Overlap parameter used for the final assembly. This is the only argument
# to the script
OL=55
# The number of threads to use
CPU=8
# To save memory, we index $D reads at a time then merge the indices together
D=4000000
# Correction k-mer value
CK=33
# The minimum k-mer coverage for the filter step. Each 27-mer
# in the reads must be seen at least this many times
COV_FILTER=2
# Overlap parameter used for FM-merge. This value must be no greater than the minimum
# overlap value you wish to try for the assembly step.
MOL=45
# Parameter for the small repeat resolution algorithm
R=10
# The number of pairs required to link two contigs into a scaffold
MIN_PAIRS=1
# The minimum length of contigs to include in a scaffold
MIN_LENGTH=100
# Distance estimate tolerance when resolving scaffold sequences
SCAFFOLD_TOLERANCE=1
# Turn off collapsing bubbles around indels
MAX_GAP_DIFF=0.05

# First, preprocess the data to remove ambiguous basecalls
sga preprocess --phred64 --pe-mode 1 -o ./preprocessed.fastq $IN1 $IN2
#
# Error correction
#
# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
sga index -a ropebwt -t $CPU --no-reverse preprocessed.fastq
# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
sga correct -k $CK --discard --learn -t $CPU -o reads.ec.k$CK.fastq preprocessed.fastq

#
# Contig assembly	
#
# Index the corrected data.
sga index -a ropebwt -t $CPU reads.ec.k$CK.fastq

# Remove exact-match duplicates and reads with low-frequency k-mers
sga filter -x $COV_FILTER -t $CPU --homopolymer-check --low-complexity-check reads.ec.k$CK.fastq

# Merge simple, unbranched chains of vertices
sga fm-merge -m $MOL -t $CPU -o merged.k$CK.fa reads.ec.k$CK.filter.pass.fa

# Build an index of the merged sequences
sga index -d 1000000 -t $CPU merged.k$CK.fa

# Remove any substrings that were generated from the merge process
sga rmdup -t $CPU merged.k$CK.fa

# Compute the structure of the string graph
sga overlap -m $MOL -t $CPU reads.ec.k$CK.filter.pass.fa

# Perform the contig assembly without bubble popping
sga assemble -m $OL -g MAX_GAP_DIFF -r $R -o assemble.m$OL reads.ec.k$CK.filter.pass.asqg.gz


#----------------ASSEMBLY STEP END----------------------------
#----------------POST-ASSEMBLY STEP----------------------------
CTGS=assemble.m$OL-contigs.fa
GRAPH=assemble.m$OL-graph.asqg.gz
sga-align --name lib.pe -t 8 $CTGS $IN1 $IN2
sga-bam2de.pl -n $MIN_PAIRS --prefix libPE lib.pe.bam
DistanceEst -s 100 --mind -99 -n 1 -k 45 -j 1 -o libPE.de libPE.hist -l 45 libPE.diffcontigs.sorted.bam
sga-astat.py -m $MIN_LENGTH lib.pe.refsort.bam > libPE.astat
sga scaffold -m $MIN_LENGTH --pe libPE.de -a libPE.astat -o scaffolds.n$MIN_PAIRS.scaf $CTGS
sga scaffold2fasta -m $MIN_LENGTH -a $GRAPH -o scaffolds.n$MIN_PAIRS.fa -d $SCAFFOLD_TOLERANCE --use-overlap --write-unplaced scaffolds.n$MIN_PAIRS.scaf
#----------------POST-ASSEMBLY STEP END----------------------------

#----------------RUN MAGNOLIA STEP----------------------------
#A bit tricky. First we need to align reads to contigs and transform them to .ace format
KAligner -m -j 8 --seq -k 61 preprocessed.fastq $CTGS > KAlign.out
perl $PATHTOABYSSTOACE/abyss2ace preprocessed.fastq $CTGS KAlign.out >abyss.ace
#This script was altered to put all reads into one group, without regexp settings
#Also additional corrections were done (output for last contig was always missing)
python ./tools/ace2magnolia.py -a abyss.ace -r "0;C;c" -t abyss -o counts.txt

mv merged* ./data/
mv reads* ./data/
mv preprocessed* ./data/
mv dotGraph.dot ./data/
mv error* ./data/