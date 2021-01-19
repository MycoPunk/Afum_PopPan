#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2 --mem 2gb

module load bcftools
if [ -f config.txt ]; then
	source config.txt
else
	echo "need a config.txt"
fi

if [ -z $FINALVCF ]; then
	echo "Need to define FINALVCF in config.txt"
	exit
fi
for TYPE in SNP INDEL
do
IN=$FINALVCF/$PREFIX.$TYPE.combined_selected.vcf.gz
STEP1=$FINALVCF/$PREFIX.$TYPE.no_missing_no_fixed.vcf.gz
STEP2=$FINALVCF/$PREFIX.$TYPE.no_missing_no_fixed.tsv
bcftools view -e 'GT="." || AF=1 || AF=0' -Oz -o $STEP1 $IN
bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT\t%AD]\n' $STEP1 > $STEP2
done

