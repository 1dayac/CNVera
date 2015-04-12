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

