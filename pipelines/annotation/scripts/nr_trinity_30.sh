#!/usr/bin/bash
#SBATCH -p short -N 1 -n 32 --mem 32gb --out logs/trinity_cdhit_cluster.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
module load cd-hit
cd-hit-est -c 0.95 -n 10 -d 0 \
-i lib/informant/Trinity_30strain.fasta -o lib/informant/Trinity_30strain.nr \
-M 32000 -T $CPU
