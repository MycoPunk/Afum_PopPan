#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 16gb --out structure.%a.log

module unload miniconda2
module load faststructure

module list

#the numbers of K to investigate is equal to the number of array jobs requested
K=${SLURM_ARRAY_TASK_ID}

#set N_REP equal to the number of iterations to run for each value of K. 
for N_REP in {1..30}; do \
mkdir rep_${N_REP}
structure.py -K $K --input=Afum_260 --output=rep_${N_REP}/Afum_260
done
