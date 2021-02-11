#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/run_BLAST.log

module load ncbi-blast/2.9.0+

mkdir results_fastas
cd genomes

##run BLAST with seq output

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

#run blast and output sequence hits
for i in *.rename.fa; do
blastn -query ../fastas/MAT1-2_Afu3g06170_nt  -db $i -evalue .0001 -outfmt "6 sseqid sseq" -word_size 10 -out $i.MAT1-2.forfasta.out
done

for i in *.rename.fa; do
blastn -query ../fastas/MAT1-1_Afum_CDSregion_Paoletti_nt  -db $i -evalue .0001 -outfmt "6 sseqid sseq" -word_size 10 -out $i.MAT1-1.forfasta.out
done

mv *MAT1-2.forfasta.out ../results_fastas
mv *MAT1-1.forfasta.out ../results_fastas

##format output files to fasta format
cd ../results_fastas
#add carrot
for i in *.forfasta.out; do
awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $i > $i.carrot.fa
done

#add newline
for i in *.scaffolds.fa.rename.fa.MAT1-1.forfasta.out.carrot.fa; do
name=$(basename "$i" .scaffolds.fa.rename.fa.MAT1-1.forfasta.out.carrot.fa)
sed 's/\t/\n/g' $i > $name.MAT1-1.temp
done

for i in *.scaffolds.fa.rename.fa.MAT1-2.forfasta.out.carrot.fa; do
name=$(basename "$i" .scaffolds.fa.rename.fa.MAT1-2.forfasta.out.carrot.fa)
sed 's/\t/\n/g' $i > $name.MAT1-2.temp
done

#add counter
for i in *.temp; do
name=$(basename "$i" .temp)
awk '/^>/{$0=$0"-"(++i)}1' $i > $name.fasta
done

mkdir MAT1-1_hits  
mkdir MAT1-2_hits   

mv *MAT1-1.fasta MAT1-1_hits
mv *MAT1-2.fasta MAT1-1_hits

#clean up
rm *carrot.fa 
rm *.forfasta.out
rm *.temp
#rm *fa.out
