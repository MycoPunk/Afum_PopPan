#!/bin/bash
#SBATCH --job-name=remove_nas
#SBATCH --output=logs/remove_nas.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

module load VCFtools

vcftools --gzvcf Afum_global_v3.LofgrenPanGenome.SNP.combined_selected.vcf.gz --max-missing 1 --recode --out Afum_global_v3.LofgrenPanGenome.SNP.combined_selected.vcf.gz

mv Afum_global_v3.LofgrenPanGenome.SNP.combined_selected.vcf.gz.recode.vcf Afum_global_v3.LofgrenPanGenome.SNP.combined_selected_NO_INV.vcf
