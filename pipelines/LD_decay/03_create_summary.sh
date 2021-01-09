#!/usr/bin/bash
#SBATCH --mem=64G -p batch --nodes 1 --ntasks 2 --out summary.log

#Clade1
for N in {1..50}; do \
cat Clade_1_samp${N}_LD_out.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > Clade_1_samp${N}_LD_out.ld.summary
done

#Clade2
for N in {1..50}; do \
cat Clade_2_samp${N}_LD_out.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > Clade_2_samp${N}_LD_out.ld.summary
done

#Clade3
cat Clade_3_LD_out.ld | sed 1,1d | awk -F " " 'function abs(v) {return v < 0 ? -v : v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' | sort -k1,1n > Clade_3_LD_out.ld.summary

#cleanup
rm *ld
