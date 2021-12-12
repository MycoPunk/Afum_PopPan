#!/bin/bash
#SBATCH --ntasks 24 --nodes 1 --mem 96G -p intel,batch
#SBATCH --time 8:00:00 --out logs/iprscan.%a.log -a 1-31
#--mail-type=END # notifications for job done & fail
#--mail-user=jasonst@ucr.edu # send-to address

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
module unload miniconda2
module load miniconda3
module load funannotate
module load iprscan/5.51-85.0
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
	N=1
    fi
fi
pushd Orthologs.longest
PREFIX=Afpan
INFILE=${PREFIX}.${N}
OUT=${INFILE}.iprout.tsv
if [ ! -s $OUT ]; then
    interproscan.sh --goterms --iprlookup -o $OUT -pa -i $INFILE -f TSV --cpu $CPU
fi
