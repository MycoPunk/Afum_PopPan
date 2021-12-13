#!/usr/bin/bash -l
#SBATCH -p short

module load miniconda3
source activate GEN220

./scripts/get_longest_protein_OrthoGroup.py

mkdir -p Orthologs.longest
bp_dbsplit.pl --prefix Orthologs.longest/Afpan --size 500 -i Orthologs.Longest.fa
