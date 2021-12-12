#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 16 --nodes 1 --mem 24G --out logs/predict.%a.log
module unload perl
module unload python
module unload miniconda2
module unload anaconda3
module load miniconda2
module load funannotate/1.8.0
source activate funannotate-1.8
which diamond
which perl
#/opt/linux/centos/7.x/x86_64/pkgs/genemarkESET/4.59_lic/gmes_petap.pl 
#export AUGUSTUS_CONFIG_PATH=/bigdata/stajichlab/shared/pkg/augustus/3.3/config
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
GMFOLDER=`dirname $(which gmhmme3)`
#genemark key is needed
if [ ! -f ~/.gm_key ]; then
	ln -s $GMFOLDER/.gm_key ~/.gm_key
fi

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
species="Aspergillus fumigatus"
IFS=,  # this species the column split character in $SAMPFILE
BUSCO=$(realpath eurotiales_odb10)
cat $SAMPFILE | sed -n ${N}p | while read BASE LOCUS
do
	name=$BASE
 	SEED_SPECIES="aspergillus_fumigatus"
	if [ ! -f $INDIR/$name.masked.fasta ]; then
		echo "No genome for $INDIR/$name.masked.fasta yet - run 00_mash.sh $N"
		exit
	fi
	PEP=$(realpath lib/informant.aa)
	OUT=$(realpath $OUTDIR/$name)
	INP=$(realpath $INDIR/$name.masked.fasta)
 	mkdir /scratch/$name.predict.$$
 	pushd /scratch/$name.predict.$$
    	funannotate predict --cpus $CPU --keep_no_stops --SeqCenter UCR --busco_db $BUSCO --strain "$BASE" \
      -i $INP --name $LOCUS --protein_evidence $PEP \
      -s "$species"  -o $OUT --busco_seed_species $SEED_SPECIES
	popd
 	rmdir /scratch/$name.predict.$$
done
