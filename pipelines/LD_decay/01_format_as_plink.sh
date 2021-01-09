#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out reformat.log

#make the bed files that structure needs using plink
module load plink

#Clade1
for N in {1..50}; do \
plink --vcf vcf/Clade_1_samp_${N}.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out Clade_1_samp${N}
done

#Clade2
for N in {1..50}; do \
plink --vcf vcf/Clade_2_samp_${N}.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out Clade_2_samp${N}
done

#Clade3
plink --vcf vcf/Clade_3.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out Clade_3
