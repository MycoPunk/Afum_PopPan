#!/usr/bin/bash
#SBATCH --mem=500G -p highmem --nodes 1 --ntasks 2 --out logs/average_summaries_all.log

module load R

Rscript pipeline/05_average_summaries_all.R
