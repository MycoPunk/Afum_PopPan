#!/bin/bash
#SBATCH -p short logs/find_missing_masked.log

CPU=1

INDIR=genomes
OUTDIR=genomes
SAMPFILE=samples.csv
IFS=,
N=1
cat $SAMPFILE | while read NAME PREF
do
 name=$NAME
 if [ ! -f $INDIR/${name}.scaffolds.fasta ]; then
#    echo -e "\tCannot find $INDIR/$name.scaffolds.fasta in $INDIR - may not have been run yet ($N)\n"
	sleep 0
 elif [ ! -f $OUTDIR/${name}.masked.fasta ]; then
#	 echo "need to run $name ($N)"
	 echo $N
 fi
 N=$(expr $N + 1)
done | perl -p -e 's/\n/,/' | perl -p -e 's/,$//'

#N=1

#m=$(cat $SAMPFILE | while read NAME PREFIX
#do
# name=$PREFIX
# if [[ -f $INDIR/${name}.sorted.fasta && ! -f $OUTDIR/${name}.masked.fasta ]]; then
#         echo $N
# fi
# N=$(expr $N + 1)
#done | perl -p -e 's/\n/,/' | perl -p -e 's/,$//')

echo "sbatch --array=$m pipeline/00_mask.sh"
