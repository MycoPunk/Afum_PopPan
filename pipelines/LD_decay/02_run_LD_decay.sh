#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/LD.log

module load plink

CPUS=$SLURM_TASKS_PER_NODE
if [ ! $CPUS ]; then
 CPUS=$PBS_NP
 if [ ! $CPUS ]; then
  CPUS=2
 fi
fi

mkdir LD_out

#Clade1
for N in {1..20}; do \
plink --bfile "plink_files/Clade_1_samp${N}" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 500 --threads $CPUS --out "LD_out/Clade_1_samp${N}_LD_out"
done

#Clade2
for N in {1..20}; do \
plink --bfile "plink_files/Clade_2_samp${N}" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 500 --threads $CPUS --out "LD_out/Clade_2_samp${N}_LD_out"
done

#Clade3
for N in {1..20}; do \
plink --bfile "plink_files/Clade_3_samp${N}" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 500 --threads $CPUS --out "LD_out/Clade_3_samp${N}_LD_out"
done

#All
for N in {1..20}; do \
plink --bfile "plink_files/n_12_samp${N}" --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 99999 --ld-window-kb 500 --threads $CPUS --out "LD_out/n_12_samp${N}_LD_out"
done

#clean up 
rm LD_out/*LD_out.log
rm LD_out/*LD_out.nosex
