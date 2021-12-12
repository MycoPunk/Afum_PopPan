#!/bin/bash
#SBATCH --nodes 1 --ntasks 16 --mem 256gb -J prepafumAAFTF --out logs/AAFTF_readprep.%a.log -p highmem --time 24:00:00

hostname
MEM=256
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

module load AAFTF/git-live

OUTDIR=input
SAMPLEFILE=samples.dat
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
PHYLUM=$(sed -n ${N}p $SAMPLEFILE | cut -f2)

if [ -z $CPU ]; then
    CPU=1
fi

WORKDIR=working_AAFTF
LEFTTRIM=$WORKDIR/${BASE}_1P.fastq.gz
RIGHTTRIM=$WORKDIR/${BASE}_2P.fastq.gz

LEFT=$WORKDIR/${BASE}_filtered_1.fastq.gz
RIGHT=$WORKDIR/${BASE}_filtered_2.fastq.gz

mkdir -p $WORKDIR

echo "$BASE"
if [ ! -f $LEFT ]; then
	echo "$OUTDIR/${BASE}_R1.fq.gz $OUTDIR/${BASE}_R2.fq.gz"
	if [ ! -f $LEFTTRIM ]; then
	    AAFTF trim --method bbduk --memory $MEM --left $OUTDIR/${BASE}_R1.fq.gz --right $OUTDIR/${BASE}_R2.fq.gz -c $CPU -o $WORKDIR/${BASE}
	fi
	echo "$LEFTTRIM $RIGHTTRIM"
	AAFTF filter -c $CPU --memory $MEM -o $WORKDIR/${BASE} --left $LEFTTRIM --right $RIGHTTRIM --aligner bbduk 
	echo "$LEFT $RIGHT"
	if [ -f $LEFT ]; then
		unlink $LEFTTRIM
		unlink $RIGHTTRIM
	fi
fi
