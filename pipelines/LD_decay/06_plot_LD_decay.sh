#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/make_plots.log

module load R

Rscript pipeline/06_plot_LD_decay.R
