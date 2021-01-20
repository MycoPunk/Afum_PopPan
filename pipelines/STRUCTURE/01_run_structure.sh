#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 16gb --out structure.%a.log

module load faststructure

module list
module unload python
module unload miniconda3
module load miniconda2

#the numbers of K to investigate is equal to the number of array jobs requested (I ran K = 1-15)
K=${SLURM_ARRAY_TASK_ID}

#set N_REP equal to the number of iterations to run for each value of K.
for N_REP in {1..30}; do \
mkdir rep_${N_REP}
structure.py -K $K --input=Afum_262 --output=rep_${N_REP}/Afum_262
done
