# !/usr/bin/env python
# coding=utf-8
# Writer:         chaijingchao
# Program Date:   2016.07.16

from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description="stat the length of fasta")
parser.add_argument('--in',help='fasta file',required=True)
parser.add_argument('--out',help='sequence length file',required=True)
parser.add_argument('--cDNA',help='cDNA size',required=True)
argv=vars(parser.parse_args())
input = open(argv["in"])
output = open(argv["out"],"w")
cDNA = argv["cDNA"]
for seq_record in SeqIO.parse(input,"fasta"):
     gene_id = seq_record.id
     seq = str(seq_record.seq).strip()
     gene_length = len(seq)
     output.write(cDNA + "\t" + gene_id + "\t" + str(gene_length) + "\n")
output.close()
