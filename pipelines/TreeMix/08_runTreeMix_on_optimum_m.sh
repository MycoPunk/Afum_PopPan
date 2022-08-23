#!/bin/bash
#SBATCH --job-name=run_TreeMix_optimum
#SBATCH --output=logs/run_TreeMix.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=4:00:00
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --partition=common

module load TreeMix/1.13

#run treemix
treemix -i Afum_global_v3.LofgrenPanGenome.SNP.combined_selected_NO_INV.vcf.LDpruned.treemix.frq.gz -root outgroup -m 1

#pack them up
mkdir TreeMix_results_optimum
mv TreeMix* TreeMix_results_optimum/
