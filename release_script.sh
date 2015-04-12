#First you need to do is simulate reads from reference
#For this purpose we used art - http://www.niehs.nih.gov/research/resources/software/biostatistics/art/  
#Please export art_illumina to paths
cd data
art_illumina -i DH10B-K12.fasta -p -f 25 -l 100 -m 300 -s 50 -qs 50 -qs2 50 -o R
#Then for read correction start Quake
echo -e "R1.fq\nR2.fq" >readnames.txt
python quake.py -f readnames.txt -k 17

cd ..