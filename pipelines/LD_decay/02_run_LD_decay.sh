#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out LD.log

#Clade1
for N in {1..50}; do \
plink --bfile "Clade_1_samp${N}" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "Clade_1_samp${N}_LD_out"
done

#Clade2
for N in {1..50}; do \
plink --bfile "Clade_2_samp${N}" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "Clade_2_samp${N}_LD_out"
done

#Clade3
plink --bfile "Clade_3" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 999999 --ld-window-kb 8000 --out "Clade_3_LD_out"

#clean up 
rm *LD_out.log
rm *LD_out.nosex
