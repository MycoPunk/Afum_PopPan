#!/bin/bash
#SBATCH -p batch --time 2-0:00:00 --ntasks 8 --nodes 1 --mem 24G --out logs/mask.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
OUTDIR=RepeatMasker_run

LIBRARY=$(realpath lib/Afum95_Fungi_repeats.lib)
SAMPFILE=samples.csv
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

IFS=,
SPECIES="Aspergillus fumigatus"
cat $SAMPFILE | sed -n ${N}p | while read BASE LOCUS
do
  name=$BASE
 if [ ! -f $INDIR/${name}.scaffolds.fasta ]; then
     echo "Cannot find $name in $INDIR - may not have been run yet"
     exit
 fi

 if [ ! -s $INDIR/${name}.masked.fasta ]; then
	 module unload miniconda2
module unload anaconda3
module load miniconda3
module unload perl
module unload python
module load funannotate/1.8.0
source activate funannotate-1.8
module switch ncbi-rmblast/2.9.0-p2

#  funannotate mask --cpus $CPU -i ../$INDIR/${name}.scaffolds.fasta -o ../$OUTDIR/${name}.masked.fasta -l $LIBRARY --method repeatmasker
  mkdir -p $OUTDIR/${name}
  RepeatMasker -e ncbi -xsmall -s -pa $CPU -lib $LIBRARY -dir $OUTDIR/${name} -gff $INDIR/${name}.scaffolds.fasta 
  rsync -a $OUTDIR/${name}/${name}.scaffolds.fasta.masked $INDIR/${name}.masked.fasta
else
    echo "Skipping ${name} as masked already"
fi

done
