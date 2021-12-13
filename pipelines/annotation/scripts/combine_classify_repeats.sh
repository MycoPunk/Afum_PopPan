#!/usr/bin/bash
#SBATCH -p short -N 1 -n 24 --mem 32gb --out logs/combine_repeats.log

CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi
CPU=$CPUS

if [ ! -s lib/repeats.all ]; then
 for n in $(ls repeat_library/*-library.fasta); do m=$(basename $n .repeatmodeler-library.fasta); perl -p -e "s/>/>$m\_/" $n; done > lib/repeats.all
fi

module load cd-hit
module load diamond
module load ncbi-blast/2.9.0+
module load miniconda3
source activate bcbio
DB=/srv/projects/db/Swissprot/2020_08/uniprot_sprot
 for n in 99 97 95 90
 do
	 if [ ! -s lib/repeats.nr_$n ]; then
 		cd-hit -T $CPU -i lib/repeats.all -o lib/repeats.nr_$n -c 0.${n}
	fi
	if [ ! -s lib/repeats.nr_${n}.diamond_blastx.swissprot.tab ]; then
	   diamond blastx -q lib/repeats.nr_${n} --db $DB -o lib/repeats.nr_${n}.diamond_blastx.swissprot.tab --outfmt 6
	fi
	if [ ! -s lib/repeats.nr_${n}.BLASTX.swissprot.tab ]; then
    #-sorthits 4 -sorthsps 2
 	  blastx -query lib/repeats.nr_${n} -db $DB  -outfmt 6 -out lib/repeats.nr_${n}.BLASTX.swissprot.tab -evalue 1e-3  -max_intron_length 1000 -max_target_seqs 10 -num_threads $CPU
	fi
  perl scripts/classify_repeats.pl -p $n -hits lib/repeats.nr_${n}.BLASTX.swissprot.tab > lib/repeats.nr_${n}.reclassify.lib
 done
