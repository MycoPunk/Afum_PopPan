#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out logs/summary.log

#Clade1
for N in {1..50}; do \
cat LD_out/Clade_1_samp${N}_LD_out.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > LD_out/Clade_1_samp${N}_LD_out.ld.summary
done

#Clade2
for N in {1..50}; do \
cat LD_out/Clade_2_samp${N}_LD_out.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > LD_out/Clade_2_samp${N}_LD_out.ld.summary
done

#Clade3
for N in {1..50}; do \
cat LD_out/Clade_3_samp${N}_LD_out.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > LD_out/Clade_3_samp${N}_LD_out.ld.summary
done

#cleanup
rm LD_out/*_out.ld
