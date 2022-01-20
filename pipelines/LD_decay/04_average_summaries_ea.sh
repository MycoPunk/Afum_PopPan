#!/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/average_summaries_ea.log

module load R

for N in LD_out/*; do \
name=$(basename "$N")
Rscript pipeline/04_average_summaries_ea.R $name --save
done
