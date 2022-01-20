#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/reformat.log

#make the bed files that structure needs using plink
module load plink

mkdir plink_files

#Clade1
for N in {1..20}; do \
plink --vcf vcf/Clade_1_samp_${N}.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out plink_files/Clade_1_samp${N}
done

#Clade2
for N in {1..20}; do \
plink --vcf vcf/Clade_2_samp_${N}.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out plink_files/Clade_2_samp${N}
done

#Clade3
for N in {1..20}; do \
plink --vcf vcf/Clade_3_samp_${N}.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out plink_files/Clade_3_samp${N}
done

#All
for N in {1..20}; do \
plink --vcf vcf/n_12_samp_${N}.recode_noinv.vcf --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out plink_files/n_12_samp${N}
done
