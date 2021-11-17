#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 16gb --out chooseK.%a.log

module unload miniconda2
module load faststructure

for N_REP in {1..30}; do \
chooseK.py --input=rep_${N_REP}/Afum_260
done
