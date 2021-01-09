#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out rairfy.log

#load modules 
module load vcftools
module load bcftools 


####
#This script takes in three clade assignment txt files (clade assignment determined previously using DAPC): Calde_1_members.txt, #Calde_2_members.txt, Calde_3_members.txt
#they are just line seperated names 
#>cat Calde_3_members.txt | head -3
#CM2495
#CM2730
#CM2733
#The script makes random subsets of names of the smallest clade size (15 in our case), subsets theses from the main vcf file, and removes invariant sites. 


#randomly sort text files of strain names and take the first 15 lines (50 times for each clade). 
#do this 50 times
for N in {1..50}; do \
   sort -R Clade_1_members.txt | head -n 15 > C1_sample_${N}.txt
done

for N in {1..50}; do \
   sort -R Clade_2_members.txt | head -n 15 > C2_sample_${N}.txt
done

#rarify the vcf files from the randomly sampled strain name lists created above
for N in {1..50}; do \
vcftools --keep C1_sample_${N}.txt --gzvcf vcf/pop_for_pan2_.SNP.combined_selected.NO_TEs.vcf.gz --recode --out vcf/Clade_1_samp_${N}
done

#rarify the vcf files from the randomly sampled name lists
for N in {1..50}; do \
vcftools --keep C2_sample_${N}.txt --gzvcf vcf/pop_for_pan2_.SNP.combined_selected.NO_TEs.vcf.gz --recode --out vcf/Clade_2_samp_${N}
done

#subset off of clade 3 (this n=15 group is what we're rairfying to).
vcftools --keep Clade_3_members.txt --gzvcf vcf/pop_for_pan2_.SNP.combined_selected.NO_TEs.vcf.gz --recode --out vcf/Clade_3

#remove invariant sites using bcftools
#Clade 1
for N in {1..50}; do \
bcftools view -c 1 vcf/Clade_1_samp_${N}.recode.vcf -o vcf/Clade_1_samp_${N}.recode_noinv.vcf
done

#Clade 2
for N in {1..50}; do \
bcftools view -c 1 vcf/Clade_2_samp_${N}.recode.vcf -o vcf/Clade_2_samp_${N}.recode_noinv.vcf
done

#Clade 3
bcftools view -c 1 vcf/Clade_3.recode.vcf -o vcf/Clade_3.recode_noinv.vcf
