#!/usr/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16 --mem 16gb
#SBATCH --output=logs/annotfunc.%a.log
#SBATCH --time=2-0:00:00
#SBATCH -p intel -J annotfunc

module unload perl
module unload python
module unload miniconda2
module unload anaconda3
module load miniconda2
module load funannotate/1.8.0
source activate funannotate-1.8
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
module load phobius
CPUS=$SLURM_CPUS_ON_NODE
OUTDIR=annotate
INDIR=genomes
SAMPFILE=samples.csv
BUSCO=eurotiomycetes_odb10
species="Aspergillus fumigatus"
if [ -z $CPUS ]; then
 CPUS=1
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=`wc -l $SAMPFILE | awk '{print $1}'`

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi
IFS=,
SPECIES="Aspergillus_fumigatus"
cat $SAMPFILE | sed -n ${N}p | while read PREFIX LOCUS
do
	name=$PREFIX
	Strain=$PREFIX
	MOREFEATURE=""
	TEMPLATE=$(realpath lib/pangenome.sbt)
	if [ ! -f $TEMPLATE ]; then
		echo "NO TEMPLATE for $name"
		exit
	fi
	ANTISMASHRESULT=$OUTDIR/$name/annotate_misc/antiSMASH.results.gbk
	echo "$name $species"
	if [[ ! -f $ANTISMASHRESULT && -d $OUTDIR/$name/antismash_local ]]; then
		ANTISMASH=$OUTDIR/$name/antismash_local/${SPECIES}_$name.gbk
		if [ ! -f $ANTISMASH ]; then
			echo "CANNOT FIND $ANTISMASH in $OUTDIR/$name/antismash_local"
		else
			rsync -a $ANTISMASH $ANTISMASHRESULT
		fi
	fi
	# need to add detect for antismash and then add that
	funannotate annotate --sbt $TEMPLATE --busco_db $BUSCO -i $OUTDIR/$name --species "$species" --strain "$Strain" --cpus $CPUS $MOREFEATURE $EXTRAANNOT
done
