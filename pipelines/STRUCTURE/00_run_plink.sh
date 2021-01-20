
#!/usr/bin/bash
#SBATCH --mem 2gb

#make the bed files that structure needs using plink

module load plink
plink --vcf vcf/262_strains.selected.SNP.NO_TEs.vcf.gz --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed --out Afum_262
