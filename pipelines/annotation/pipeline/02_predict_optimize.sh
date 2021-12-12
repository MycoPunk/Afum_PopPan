#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 16 --nodes 1 --mem 24G --out logs/predict.%a.log
# Define program name
PROGNAME=$(basename $0)

# Load software
module unload miniconda2
module unload anaconda3
module load miniconda3
module unload perl
module unload python
module load funannotate/1.8.0
source activate funannotate-1.8

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
OUTDIR=annotate

SAMPFILE=samples.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)
#export GENEMARK_PATH=/opt/genemark/gm_et_linux_64

export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
# make genemark key link required to run it
if [ ! -f ~/.gm_key ]; then
	module  load    genemarkESET/4.59_lic
	GMFOLDER=`dirname $(which gmhmme3)`
        ln -s $GMFOLDER/.gm_key ~/.gm_key
fi

IFS=,
SPECIES="Aspergillus fumigatus"
RNASEQSET=PRJNA376829
BUSCO=eurotiomycetes_odb9
cat $SAMPFILE | sed -n ${N}p | while read BASE LOCUS
do
    STRAIN=$BASE
    MASKED=$(realpath $INDIR/$BASE.masked.fasta)
    if [ ! -f $MASKED ]; then
	echo "Cannot find $BASE.masked.fasta in $INDIR - may not have been run yet"
	exit
    fi
    SEED_SPECIES="aspergillus_fumigatus"
    echo "looking for $MASKED to run"

    mkdir $BASE.predict.$$
    pushd $BASE.predict.$$
    funannotate predict --cpus $CPU --keep_no_stops --SeqCenter UCR --busco_db $BUSCO --optimize_augustus \
	--strain $STRAIN --min_training_models 100 --AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH \
	-i ../$INDIR/$BASE.masked.fasta --name $LOCUS --protein_evidence ../lib/informant.aa \
	-s "$SPECIES"  -o ../$OUTDIR/$BASE --busco_seed_species $SEED_SPECIES
    popd
    rmdir $BASE.predict.$$
done
