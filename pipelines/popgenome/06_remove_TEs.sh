#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 2gb

###use bedtools subtract to remove known TE regions

cd vcf/

module load bedtools
module load tabix

#SNPS
bedtools subtract -a *SNP.combined_selected.vcf.gz -b ../genome/FungiDB-39_AfumigatusAf293_Genome.RM.tab -header > 262_strains.selected.SNP.NO_TEs.vcf 
#INDELS
bedtools subtract -a *SNP.combined_selected.vcf.gz -b ../genome/FungiDB-39_AfumigatusAf293_Genome.RM.tab -header > 262_strains.selected.INDEL.NO_TEs.vcf 


#compress
bgzip 262_strains.selected.SNP.NO_TEs.vcf
bgzip 262_strains.selected.INDEL.NO_TEs.vcf

#index:
tabix -p vcf 262_strains.selected.SNP.NO_TEs.vcf.gz
tabix -p vcf 262_strains.selected.INDEL.NO_TEs.vcf.gz
