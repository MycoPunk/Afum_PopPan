#!/bin/bash
#SBATCH --nodes 1 --ntasks 4 --mem 16G --time 24:00:00 --out logs/busco.%a.log -J busco

module load busco
module unload ncbi-blast
module load ncbi-blast/2.2.31+

# for augustus training
# set to a local dir to avoid permission issues and pollution in global
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)

CPU=${SLURM_CPUS_ON_NODE}
N=${SLURM_ARRAY_TASK_ID}
if [ ! $CPU ]; then
     CPU=2
fi

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
if [ -z ${SLURM_ARRAY_JOB_ID} ]; then
	SLURM_ARRAY_JOB_ID=$$
fi
GENOMEFOLDER=genomes_scaffolded
EXT=scaffolds.fasta
LINEAGE=/srv/projects/db/BUSCO/v9/eurotiomycetes_odb9
OUTFOLDER=BUSCO
mkdir -p $OUTFOLDER
TEMP=/scratch/${SLURM_ARRAY_JOB_ID}_${N}
mkdir -p $TEMP
SAMPLEFILE=samples.dat
STRAIN=$(sed -n ${N}p $SAMPLEFILE | cut -f2)
SEED_SPECIES=aspergillus_fumigatus
GENOMEFILE=$(realpath $GENOMEFOLDER/${STRAIN}.${EXT})
LINEAGE=$(realpath $LINEAGE)
NAME=$(echo "$STRAIN" | perl -p -e 's/[\)\(]//g')
if [ -d "$OUTFOLDER/run_${NAME}" ];  then
    echo "Already have run $NAME in folder busco - do you need to delete it to rerun?"
    exit
else
    pushd $OUTFOLDER
    run_BUSCO.py -i $GENOMEFILE -l $LINEAGE -o $NAME -m geno --cpu $CPU --tmp $TEMP -sp $SEED_SPECIES
    popd
fi

rm -rf $TEMP
