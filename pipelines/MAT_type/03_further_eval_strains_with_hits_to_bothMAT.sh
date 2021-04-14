#!/usr/bin/bash
#SBATCH --mem 10G --ntasks 20 --nodes 1  --time 10:00:00 

#note - first symlink the R1 and R2 raw reads for the strains in question into raw_fastas folder

module load bowtie2/2.3.4.1 
module load samtools/1.11
module load bcftools
module load picard/2.18.3
module load gatk/4

cd raw_fastas/

#index the ‘references’ uses the format : bowtie2-build -f [reference file] [output name]
bowtie2-build -f ../fastas/MAT1-1_AY898661.1_whole_gene.nt MAT1-1_REF
bowtie2-build -f ../fastas/MAT1-2_Afu3g06170_whole_gene.nt MAT1-2_REF

##make SAM alignment files
#MAT1
for i in $(ls *.fq.gz | rev | cut -c 10- | rev | uniq); do
	bowtie2 --very-sensitive-local -p 20 -x MAT1-1_REF -1 $i\_R1.fq.gz -2 $i\_R2.fq.gz -S $i.MAT1.sam ;
done

#MAT2
for i in $(ls *.fq.gz | rev | cut -c 10- | rev | uniq); do
	bowtie2 --very-sensitive-local -p 20 -x MAT1-2_REF -1 $i\_R1.fq.gz -2 $i\_R2.fq.gz -S $i.MAT2.sam ;
done

##convert SAM to BAM
#MAT1
for i in *.MAT1.sam ; do
name=$(basename "$i" .sam)
	samtools view -b $i > $name.bam ;
done

#MAT1
for i in *.MAT2.sam ; do
name=$(basename "$i" .sam)
	samtools view -b $i > $name.bam ;
done

#Clean up large files
#rm*.sam

##sort BAM files
#MAT1
for i in *.MAT1.bam ; do
name=$(basename "$i" .bam)
	samtools sort -T $i.temp -o $name.sorted.bam -O bam $i ;
done

#MAT2
for i in *.MAT2.bam ; do
name=$(basename "$i" .bam)
	samtools sort -T $i.temp -o $name.sorted.bam -O bam $i ;
done

#get coverage of bam files 
for i in *.sorted.bam ; do
name=$(basename "$i" .sorted.bam)
	samtools depth -aa $i > $name.coverage ;
done


##Call variants 
#MAT1
for i in *.MAT1.sorted.bam ; do
name=$(basename "$i" .sorted.bam)
	bcftools mpileup -Ou -f ../fastas/MAT1-1_AY898661.1_whole_gene.nt $i | bcftools call -m -Oz -o $name.vcf.gz ;
done

for i in *.MAT2.sorted.bam ; do
name=$(basename "$i" .sorted.bam)
	bcftools mpileup -Ou -f ../fastas/MAT1-2_Afu3g06170_whole_gene.nt $i | bcftools call -m -Oz -o $name.vcf.gz ;
done

#clean up 
#rm *.bam

#index the vcf files
for i in *.vcf.gz ; do
tabix $i
done

#re-index the reference
#note - the ref needs to have the extension ".fasta"
cp ../fastas/MAT1-1_AY898661.1_whole_gene.nt MAT1-1_AY898661.1_whole_gene.fasta
samtools faidx MAT1-1_AY898661.1_whole_gene.fasta
cp ../fastas/MAT1-2_Afu3g06170_whole_gene.nt MAT1-2_Afu3g06170_whole_gene.fasta
samtools faidx MAT1-2_Afu3g06170_whole_gene.fasta

#create reference sequence dictionary using PicardTools
picard CreateSequenceDictionary R=MAT1-1_AY898661.1_whole_gene.fasta  O=MAT1-1_AY898661.1_whole_gene.dict
picard CreateSequenceDictionary R=MAT1-2_Afu3g06170_whole_gene.fasta  O=MAT1-2_Afu3g06170_whole_gene.dict

#print consensus fastas
for i in *.MAT1.vcf.gz ; do
cat MAT1-1_AY898661.1_whole_gene.fasta | bcftools consensus $i -a N > $i.Ns_and_variants_added.fa
done

for i in *.MAT2.vcf.gz ; do
cat MAT1-2_Afu3g06170_whole_gene.fasta | bcftools consensus $i -a N > $i.Ns_and_variants_added.fa
done
