#!/usr/bin/bash
#SBATCH -p short  --mem 16gb -N 1 -n 8 --out logs/ragtag.%a.log

module unload miniconda3
module unload miniconda3
module load anaconda3
module load minimap2
module load ragtag
conda activate ragtag

CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
if [ -z $CPU ]; then
	CPU=1
fi
ASM=genomes
OUT=scaffolds
CTGS=$(cut -f2 samples.dat | sed -n ${N}p)
CTGS=$ASM/$CTGS.sorted.fasta

NAME=$(basename $CTGS .sorted.fasta)
CTGS=$(realpath $CTGS)
REF=$(realpath genome_ref/FungiDB-39_AfumigatusAf293_Genome.fasta)
echo $REF
echo $CTGS
SCAF=$OUT/$NAME
mkdir -p $SCAF
if [ ! -s $SCAF/assembly.stats.txt ]; then
  pushd $OUT
  ragtag.py scaffold $REF $CTGS --aligner $(which minimap2) -u -o $NAME -t $CPU -r
  module load AAFTF
  AAFTF assess -i $NAME/ragtag.scaffolds.fasta -r $NAME/assembly.stats.txt
fi
