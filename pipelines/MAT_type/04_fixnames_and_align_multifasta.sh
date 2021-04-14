#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out rename_and_align.log

module load mafft/7.471

cd raw_fastas/

#attach strain names to fasta headers to keep track 
for i in *.Ns_and_variants_added.fa; do
name=$(basename "$i" .vcf.gz.Ns_and_variants_added.fa)
awk '/>/{sub(">","&"FILENAME"_from_");sub(/\.vcf.gz.Ns_and_variants_added.fa/,x)}1' $i > $name.consensus.fa
done

#move results based on MAT type
mkdir MAT1-1_hits
mkdir MAT1-2_hits

mv *MAT1.consensus.fa MAT1-1_hits/
mv *MAT2.consensus.fa MAT1-2_hits/

cd MAT1-1_hits/

#concat fastas
cat *consensus.fa > input.fasta

#run alignment
mafft input.fasta >MAT1-1_alignment.fasta

cd ../MAT1-2_hits    

#concat fastas
cat *consensus.fa > input.fasta

#run alignment
mafft input.fasta >MAT1-2_alignment.fasta
