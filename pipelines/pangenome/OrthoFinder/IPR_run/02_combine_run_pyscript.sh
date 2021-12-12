#!/usr/bin/bash
#SBATCH -p short --mem 8gb

cat Orthologs.longest/*iprout.tsv | sort > Af.pan.genome.iprout.tsv

#python3 04_combine_ipr_add_genefamily.py
