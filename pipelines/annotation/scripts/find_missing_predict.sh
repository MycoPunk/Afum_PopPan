#!/bin/bash
#SBATCH -p short logs/find_missing_masked.log

CPU=1

INDIR=annotate
SAMPFILE=samples.csv
IFS=,
N=1
m=$(cat $SAMPFILE | while read NAME PREF
do
 name=$NAME
 if [ ! -d $INDIR/${name}/predict_results ]; then
	 echo $N
 fi
 N=$(expr $N + 1)
done | perl -p -e 's/\n/,/' | perl -p -e 's/,$//')

echo "sbatch --array=$m pipeline/02_predict_optimize.sh"
