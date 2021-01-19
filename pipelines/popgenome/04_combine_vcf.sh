#!/usr/bin/bash
#SBATCH -p intel --mem 64gb -N 1 -n 4 --out logs/concat_vcf.log -p short

module load bcftools

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
if [ -f config.txt ]; then
	source config.txt
else
	echo "need a config.txt"
fi

if [ -z $FINALVCF ]; then
	echo "Need to define FINALVCF in config.txt"
	exit
fi
if [ -z $SLICEVCF ]; then
	echo "Need to define SLICEVCF in config.txt"
	exit
fi
IN=$SLICEVCF/$PREFIX
mkdir $FINALVCF
for TYPE in SNP INDEL
do
  	OUT=$FINALVCF/$PREFIX.$TYPE.combined_selected.vcf.gz
	bcftools concat -Oz -o $OUT --threads $CPU $IN.*.$TYPE.selected.vcf.gz
	tabix $OUT
done
