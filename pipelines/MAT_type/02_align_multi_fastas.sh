#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/format_fastas.log

module load mafft/7.471

cd results_fastas/MAT1-1_hits

#concat fastas
cat *fasta > input.fasta

#run alignment
mafft input.fasta >MAT1-1_alignment.fasta

cd ../MAT1-2_hits    

#concat fastas
cat *fasta > input.fasta

#run alignment
mafft input.fasta >MAT1-2_alignment.fasta
