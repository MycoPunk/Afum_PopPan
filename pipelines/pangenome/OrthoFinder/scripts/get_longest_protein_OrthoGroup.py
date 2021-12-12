#!/usr/bin/env python3

from Bio import SeqIO
import os

reffolder=""
reffoldertxt='current_orthorun.txt';
with open(reffoldertxt,'r') as ifh:
    for n in ifh:
        reffolder=n.strip()
outfile="Orthologs.Longest.fa"
longest_seqs = []
for filename in os.listdir(reffolder):
    if not filename.endswith(".fa"):
        continue
    orthogroup = filename.split(".")[0]
    filename= os.path.join(reffolder,filename)
    longest = None
    for seq_record in SeqIO.parse(filename, "fasta"):
        if not longest or len(seq_record) > len(longest):
            longest = seq_record

    if longest:
        longest.description = longest.id
        longest.id = orthogroup
        longest_seqs.append(longest)

SeqIO.write(longest_seqs,outfile,'fasta')
