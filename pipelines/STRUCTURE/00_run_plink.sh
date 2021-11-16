#!/usr/bin/bash
#SBATCH --mem 2gb

#make the bed files that structure needs using plink

module load plink
plink --vcf vcf/Pop_for_pan_260.All.SNP.combined_selected.NO.TE.vcf.gz --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out Afum_260
