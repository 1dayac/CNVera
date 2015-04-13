__author__ = 'Dmitrii'
import sys
from Bio import SeqIO

f = open(str(sys.argv[1]), 'r')
handle = open(str(sys.argv[2]), "rU")
record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

f2 = open("new_counts.txt", 'w')

for row in f.readlines():
    splitted_row = row.split()
    splitted_row[1] = str(len(record_dict[splitted_row[0]].seq))
    f2.write("\t".join(splitted_row))
    f2.write('\n')