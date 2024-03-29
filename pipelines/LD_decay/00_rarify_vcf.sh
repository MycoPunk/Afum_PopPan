#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/rairfy.log

#load modules 
module load vcftools
module load bcftools 
mkdir temp

#randomly sort text files of strain names and take the first 12 lines (20 times for each clade). 
for N in {1..20}; do \
   sort -R Clade_1_members.txt | head -n 12 > temp/C1_sample_temp${N}.txt
done

for N in {1..20}; do \
   sort -R Clade_2_members.txt | head -n 12 > temp/C2_sample_temp${N}.txt
done

for N in {1..20}; do \
   sort -R Clade_3_members.txt | head -n 12 > temp/C3_sample_temp${N}.txt
done

for N in {1..30}; do \
   sort -R all_members.txt | head -n 12 > temp/n_12_sample_temp${N}.txt
done


#subset the strains using the randomly sampled strain name lists created above, from the vcf files 
for N in {1..20}; do \
vcftools --keep temp/C1_sample_temp${N}.txt --gzvcf vcf/Pop_for_pan_260.v2.All.SNP.combined_selected.NO.TE.vcf.gz --recode --out temp/Clade_1_samp_${N}
done

for N in {1..20}; do \
vcftools --keep temp/C2_sample_temp${N}.txt --gzvcf vcf/Pop_for_pan_260.v2.All.SNP.combined_selected.NO.TE.vcf.gz --recode --out temp/Clade_2_samp_${N}
done

for N in {1..20}; do \
vcftools --keep temp/C3_sample_temp${N}.txt --gzvcf vcf/Pop_for_pan_260.v2.All.SNP.combined_selected.NO.TE.vcf.gz --recode --out temp/Clade_3_samp_${N}
done

for N in {1..20}; do \
vcftools --keep temp/n_12_sample_temp${N}.txt --gzvcf vcf/Pop_for_pan_260.v2.All.SNP.combined_selected.NO.TE.vcf.gz --recode --out temp/n_12_samp_${N}
done


#remove invariant sites using bcftools
for N in {1..20}; do \
bcftools view -c 1 temp/Clade_1_samp_${N}.recode.vcf -o vcf/Clade_1_samp_${N}.recode_noinv.vcf
done

for N in {1..20}; do \
bcftools view -c 1 temp/Clade_2_samp_${N}.recode.vcf -o vcf/Clade_2_samp_${N}.recode_noinv.vcf
done

for N in {1..20}; do \
bcftools view -c 1 temp/Clade_3_samp_${N}.recode.vcf -o vcf/Clade_3_samp_${N}.recode_noinv.vcf
done

for N in {1..20}; do \
bcftools view -c 1 temp/n_12_samp_${N}.recode.vcf -o vcf/n_12_samp_${N}.recode_noinv.vcf
done

#clean up 
rm temp/*temp*
rm temp/*recode.vcf
