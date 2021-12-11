#!/usr/bin/bash
#SBATCH -p short -N 1 -n 4 --mem 16gb

pigz -dc PRJNA376829_Bagasse_R1.fq.gz | perl -p -e 's/^(\@)(SRR\S+)\s+(\S+)/$1$3 $2/' | pigz -c > PRJNA376829_R1.fq.gz
pigz -dc PRJNA376829_Bagasse_R2.fq.gz | perl -p -e 's/^(\@)(SRR\S+)\s+(\S+)/$1$3 $2/' | pigz -c > PRJNA376829_R2.fq.gz

