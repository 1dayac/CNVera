#! /bin/bash -x

#Add structural variations to reference
./tools/SVSim/SVSim ./data/SVInput.txt ./data/DH10B-K12.fasta
#First you need to do is simulate reads from reference with SV
#For this purpose we used art - http://www.niehs.nih.gov/research/resources/software/biostatistics/art/  
#Please export art_illumina to paths
art_illumina -i ./data/DH10B-K12.sim.fasta -p -f 25 -l 100 -m 300 -s 50 -qs 50 -qs2 50 -o ./data/R
#Then for read correction start Quake
echo -e "./data/R1.fq\n./data/R2.fq" >./data/readnames.txt
count_qmers -f ./data/readnames.txt -k 17 -q 64
#as we know coverage lets chose value like 6, because quake often fails to di it by itself
correct -f ./data/readnames.txt -k 17 -c 6 -m readnames.txt.qcts -p 4
#Reads are stored in ./data/R1.cor.fq and ./data/R2.cor.fq 
rm kmers.txt error_model* r.log readnames.txt.qcts

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
sga preprocess --phred64 --pe-mode 1 -o ./data/preprocessed.fastq $IN1 $IN2
#
# Error correction
#
# Build the index that will be used for error correction
# As the error corrector does not require the reverse BWT, suppress
# construction of the reversed index
sga index -a ropebwt -t $CPU --no-reverse ./data/preprocessed.fastq

# Perform error correction with a 41-mer.
# The k-mer cutoff parameter is learned automatically
sga correct -k $CK --discard --learn -t $CPU -o ./data/reads.ec.k$CK.fastq ./data/preprocessed.fastq

#
# Contig assembly
#
# Index the corrected data.
sga index -a ropebwt -t $CPU ./data/reads.ec.k$CK.fastq

# Remove exact-match duplicates and reads with low-frequency k-mers
sga filter -x $COV_FILTER -t $CPU --homopolymer-check --low-complexity-check ./data/reads.ec.k$CK.filter.pass.fa

# Merge simple, unbranched chains of vertices
sga fm-merge -m $MOL -t $CPU -o ./data/merged.k$CK.fa ./data/reads.ec.k$CK.filter.pass.fa

# Build an index of the merged sequences
sga index -d 1000000 -t $CPU ./data/merged.k$CK.fa

# Remove any substrings that were generated from the merge process
sga rmdup -t $CPU ./data/merged.k$CK.fa

# Compute the structure of the string graph
sga overlap -m $MOL -t $CPU ./data/reads.ec.k$CK.filter.pass.fa

# Perform the contig assembly without bubble popping
sga assemble -m $OL -g MAX_GAP_DIFF -r $R -o ./data/assemble.m$OL ./data/reads.ec.k$CK.filter.pass.asqg.gz
