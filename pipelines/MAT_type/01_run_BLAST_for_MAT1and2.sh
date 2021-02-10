#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/run_BLAST.log

module load ncbi-blast/2.9.0+

cd genomes

#attach strain names to keep track 
for i in *.fa; do
name=$(basename "$i" genomes/)
awk '/>/{sub(">","&"FILENAME"_");sub(/\.fasta/,x)}1' $i > $i.rename.fa
done

#make blast db
for i in *rename.fa; do
name=`echo $i | awk -F "." '{print $1}'`
makeblastdb -dbtype nucl -in $i -title $i
done

#move up one directory
cd ..
mkdir results

#run blast and output sequence hits
for i in genomes/*.rename.fa; do
blastn -query fastas/MAT1-2_Afu3g06170_nt  -db $i -max_hsps 1 -evalue .0001 -outfmt 7 -word_size 10 -out $i.MAT1-2.out
done

for i in genomes/*.rename.fa; do
blastn -query fastas/MAT1-1_Afum_CDSregion_Paoletti_nt  -db $i -max_hsps 1 -evalue .0001 -outfmt 7 -word_size 10 -out $i.MAT1-1.out
done
#blastn -query fastas/MAT1-1_Afum_CDSregion_Paoletti_nt  -db $i -evalue .0001 -outfmt "6 sseqid sseq" -word_size 10 -out $i.MAT1-1.out

mv genomes/*.out results/ 
