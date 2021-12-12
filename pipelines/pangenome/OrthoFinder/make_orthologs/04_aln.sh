#!/usr/bin/bash
#SBATCH -p short -N 1 -n 96 --mem 64gb -C ryzen --out logs/msa_aln_mafft.log

module load mafft

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
	CPU=$SLURM_CPUS_ON_NODE
fi

make -f scripts/makefile.makealn -j $CPU
