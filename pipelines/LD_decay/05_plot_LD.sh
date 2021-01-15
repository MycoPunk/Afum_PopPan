#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/average_summaries.log

module load R

Rscript pipeline/05_plot_LD_decay_zoomed_in.R

Rscript pipeline/05_plot_LD_decay_LD50.R
