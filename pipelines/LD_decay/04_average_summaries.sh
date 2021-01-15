#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/average_summaries.log

module load R

Rscript pipeline/04_average_summaries.R
