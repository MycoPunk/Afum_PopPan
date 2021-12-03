#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/run_BLATSP.log

#this script assigns each gene in the Afu293 reference genome to OGs derrived using orthofinder. Here, we use the longest representative sequece in each orthogroup as the input sequnce to make$


module load ncbi-blast/2.12.0+

makeblastdb -in sm_input/Orthologs.Longest.fa -dbtype prot

#note, we're setting the max_target_seqs to a redic high number here and filtering results afterwards to avoid known problems with the command line version of max_target_seqs ecxluding top hit$
blastp -query sm_input/Af293_pep_short.fa  \
-db sm_input/Orthologs.Longest.fa   \
-out sm_input/OF_blastP_results.txt \
-outfmt 6 \
-max_target_seqs 500

#subset to get only the top hit - note these are already sorted in order of likelyhood (alignment length and eval)
#so we just want the lines where it's the first instance of the value in the first col (Af293 gene hits)
awk '!($1 in a){a[$1];print}' < sm_input/OF_blastP_results.txt > sm_input/OF_blastP_results_filtered.txt

#grab just the first two columns for analysis
cat sm_input/OF_blastP_results_filtered.txt | grep '^Afu*' > sm_input/temp.txt
cut -f 1,2 sm_input/temp.txt > sm_input/one_to_one_blast_hits.txt

rm sm_input/temp.txt
