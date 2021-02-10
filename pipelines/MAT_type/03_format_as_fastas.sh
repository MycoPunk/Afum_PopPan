#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/run_BLAST.log

module load ncbi-blast/2.9.0+

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
blastn -query ../fastas/MAT1-2_Afu3g06170_nt  -db $i -max_hsps 1 -evalue .0001 -outfmt "6 sseqid sseq" -word_size 10 -out $i.MAT1-2.forfasta.out
done

for i in *.rename.fa; do
blastn -query ../fastas/MAT1-1_Afum_CDSregion_Paoletti_nt  -db $i -max_hsps 1 -evalue .0001 -outfmt "6 sseqid sseq" -word_size 10 -out $i.MAT1-1.forfasta.out
done

mv *MAT1-2.forfasta.out ../results
mv *MAT1-1.forfasta.out ../results

##format output files to fasta format
cd ../results
#add carrot
for i in *.forfasta.out; do
awk '{ if ($0 ~ /_/) { printf ">"; } print $0; }' $i > $i.carrot.fa
done

#add newline
for i in *.scaffolds.fa.rename.fa.MAT1-1.forfasta.out.carrot.fa; do
name=$(basename "$i" .scaffolds.fa.rename.fa.MAT1-1.forfasta.out.carrot.fa)
sed 's/\t/\n/g' $i > $name.MAT1-1.fasta
done

for i in *.scaffolds.fa.rename.fa.MAT1-2.forfasta.out.carrot.fa; do
name=$(basename "$i" .scaffolds.fa.rename.fa.MAT1-2.forfasta.out.carrot.fa)
sed 's/\t/\n/g' $i > $name.MAT1-2.fasta
done

#add counter - note only necessary if you DON'T use -max_hsps 1 and want to make alignments for all hits
#for i in *.MAT.temp; do
#name=$(basename "$i" .MAT.temp)
#awk '/^>/{$0=$0"-"(++i)}1' $i > $name.MAT.fasta
#done

#clean up
rm *carrot.fa 
rm *.forfasta.out
#rm *MAT.temp
#rm *fa.out
