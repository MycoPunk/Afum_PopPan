#!/bin/bash
#SBATCH --nodes 1 --ntasks 4 --mem 8g --out logs/minimap.%a.log

module load minimap2

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

DB=../genome/Exophiala_dermatidis_V1.fasta
OUTDIR=strain_compare
mkdir -p $OUTDIR
ASM=asm
SAMPLEFILE=samples.dat
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
if [[ "$BASE" =~ ^\# ]]; then
 echo "skipping $BASE"
 exit
fi

SORTED=$ASM/${BASE}.sorted.fasta

minimap2 -t $CPU -cx asm5 --cs $DB $SORTED > $OUTDIR/${BASE}.minimap.paf
