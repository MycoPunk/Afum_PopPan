#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/run_BLAST.log

module load ncbi-blast/2.9.0+

cd genomes

#attach strain names to keep track 
#for i in *.fa; do
#name=$(basename "$i" genomes/)
#awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' $i > $i.rename.fa
#done

#make blast db
for i in *rename.fa; do
name=`echo $i | awk -F "." '{print $1}'`
makeblastdb -dbtype nucl -in $i -title $i
done

#move up one directory
cd ..
mkdir results

#run blast and output sequence hits
#note - removed -max_hsps 1
for i in genomes/*.rename.fa; do
blastn -query fastas/MAT1-1_AY898661.1_whole_gene.nt  -db $i -evalue .0001 -outfmt 7 -word_size 10 -out $i.MAT1-1.out
done

for i in genomes/*.rename.fa; do
blastn -query fastas/MAT1-2_Afu3g06170_whole_gene.nt  -db $i -evalue .0001 -outfmt 7 -word_size 10 -out $i.MAT1-2.out
done

mv genomes/*.out results/ 
