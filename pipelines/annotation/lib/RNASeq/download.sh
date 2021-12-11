#!/usr/bin/bash
#SBATCH -p stajichlab -N 1 -n 1 --mem 8gb --out dump.%a.log
module load sratoolkit
module load aspera
module unload perl
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
ACC=acc.txt
ACCNO=$(sed -n ${N}p $ACC | cut -f1)
ASPERA=$(which ascp)
if [ ! -f ${ACCNO}_1.fastq.gz ]; then
	prefetch --ascp-path "$ASPERA|$ASPERAKEY" $ACCNO
	fastq-dump --tmpdir /scratch --defline-seq '@$sn[_$rn]/$ri' --split-files $ACCNO
	perl -i -p -e 's/_reverse//' ${ACCNO}_2.fastq
	perl -i -p -e 's/_forward//' ${ACCNO}_1.fastq
	pigz ${ACCNO}_[12].fastq

fi
