#!/usr/bin/bash

#SBATCH --nodes 1 --ntasks 24 --mem 24G -p short -J asmCount --out logs/bbcount.%a.log --time 2:00:00
module load BBMap
hostname
MEM=24
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
SAMPLEFILE=samples.dat
BASE=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
READBASE=$(sed -n ${N}p $SAMPLEFILE | cut -f1)
ASM=genomes
OUTDIR=mapping_report
SORTED=$(realpath $ASM/${BASE}.sorted.fasta)

LEFT=$(realpath $INDIR/${READBASE}_R1.fq.gz)
RIGHT=$(realpath $INDIR/${READBASE}_R2.fq.gz)
mkdir -p $OUTDIR
if [[ ! -s $OUTDIR/${BASE}.bbmap_covstats.txt || $SORTED -nt $OUTDIR/${BASE}.bbmap_covstats.txt ]]; then
	mkdir -p N$N.$$.bbmap
	pushd N$N.$$.bbmap
	if [ ! -e $RIGHT ]; then
		bbmap.sh -Xmx${MEM}g ref=$SORTED in=$LEFT covstats=../$OUTDIR/${BASE}.bbmap_covstats.txt  statsfile=../$OUTDIR/${BASE}.bbmap_summary.txt
	else
		bbmap.sh -Xmx${MEM}g ref=$SORTED in=$LEFT in2=$RIGHT covstats=../$OUTDIR/${BASE}.bbmap_covstats.txt  statsfile=../$OUTDIR/${BASE}.bbmap_summary.txt
	fi
	popd
	rm -rf N$N.$$.bbmap
fi
