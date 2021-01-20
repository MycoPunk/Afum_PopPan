#!/usr/bin/bash
module load faststructure

for N_REP in {1..30}; do \
chooseK.py --input=rep_${N_REP}/Afum_262
done
