# !/usr/bin/env python
# coding=utf-8
# Writer:         chaijingchao
# Program Date:   2016.07.16

from Bio import SeqIO
import sys
if len(sys.argv) != 5:
        print"python iso_seq_id.py indir sample cDNA_size outdir"
        exit(0)
out1 = open(sys.argv[4] + "/high_qv_consensus_isoforms.fasta","w")
out2 = open(sys.argv[4] + "/low_qv_consensus_isoforms.fasta","w")
out3 = open(sys.argv[4] + "/id.xls","w")
dic = {}
i = 1
for seq_record in SeqIO.parse(sys.argv[1] + "/polished_high_qv_consensus_isoforms.fasta","fasta"):
        gene_id = seq_record.id.strip()
        seq = str(seq_record.seq).strip()
        out1.write(">" + sys.argv[2] + "_" + sys.argv[3] + "_" + str(i) + "\n" + seq + "\n")
        out3.write(sys.argv[2] + "_" + sys.argv[3] + "_" + str(i) + "\t" + gene_id + "\n")
        i +=1
for seq_record in SeqIO.parse(sys.argv[1] + "/polished_low_qv_consensus_isoforms.fasta","fasta"):
        gene_id = seq_record.id.strip()
        seq = str(seq_record.seq).strip()
        out2.write(">" + sys.argv[2] + "_" + sys.argv[3] + "_" + str(i) + "\n" + seq + "\n")
        out3.write(sys.argv[2] + "_" + sys.argv[3] + "_" + str(i) + "\t" + gene_id + "\n")
        i +=1
out1.close
out2.close
out3.close
