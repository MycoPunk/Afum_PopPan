#!/bin/bash

#SBATCH --job-name=prune_LD     
#SBATCH --output=logs/prune_high_LD_sites.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --partition=common

module load VCFtools
module load Plink/2.00a2LM #called w. plink2 on DCC

#read in files
file=Afum_global_v3.LofgrenPanGenome.SNP.combined_selected_NO_INV.vcf

#set threashold value
thresh=0.2

#output a map o linkage to identify which shoul be pruned
vcftools --vcf $file --plink --out $file.plink

#replace first col of 0s with 1s for formatting needs
sed -i 's/^0\t/1\t/g' ${file}.plink.map

#convert vcf to plink format
plink2 --vcf $file --set-all-var-ids @:#[b37] --double-id --allow-extra-chr --indep-pairwise 50 10 $thresh --out treemix_ldfilter

#remove the ":" for formatting
sed -i 's/:/\t/g' treemix_ldfilter.prune.in

#output pruned file
vcftools --vcf $file  --out $file.pruning --positions treemix_ldfilter.prune.in --stdout --recode | gzip > $file.LDpruned.vcf.gz
