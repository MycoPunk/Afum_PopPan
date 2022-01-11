#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/Jellyfish.log

module load jellyfish

jellyfish count -C -m 21 -s 1000000000 -t 10 genomes/Aspergillus_fumigatus_Afu_343_P_11.scaffolds.fa -o reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/IFM_61407_R1.test.fq  -o IFM_61407_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/IFM_59359_R1.test.fq  -o IFM_59359_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/DMC2_AF100-1_3_R1.test.fq  -o DMC2_AF100-1_3_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/DMC2_AF100-1_18_R1.test.fq  -o DMC2_AF100-1_18_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/08-36-03-25_R1.test.fq  -o 08-36-03-25_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/NCPF-7816_R1.test.fq  -o NCPF-7816_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/Afu_343-P-11_R1.test.fq  -o Afu_343-P-11_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/B7586_CDC-30_R1.test.fq  -o B7586_CDC-30_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/AF100-12_5_R1.test.fq  -o AF100-12_5_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/SF2S9_R1.fastq  -o SF2S9_reads.jf
jellyfish count -C -m 21 -s 1000000000 -t 10 raw_fastas/AF100-12_7G_R1.fastq  -o AF100-12_7G_reads.jf



#jellyfish histo -t 10 reads.jf > reads.histo
jellyfish histo -t 10 IFM_61407_reads.jf > IFM_61407_reads.histo
jellyfish histo -t 10 IFM_59359_reads.jf > IFM_59359_reads.histo
jellyfish histo -t 10 DMC2_AF100-1_3_reads.jf > DMC2_AF100-1_3_reads.histo
jellyfish histo -t 10 DMC2_AF100-1_18_reads.jf > DMC2_AF100-1_18_reads.histo
jellyfish histo -t 10 08-36-03-25_reads.jf > 08-36-03-25_reads.histo
jellyfish histo -t 10 NCPF-7816_reads.jf > NCPF-7816_reads.histo
jellyfish histo -t 10 Afu_343-P-11_reads.jf > Afu_343-P-11_reads.histo
jellyfish histo -t 10 B7586_CDC-30_reads.jf > B7586_CDC-30_reads.histo
jellyfish histo -t 10 AF100-12_5_reads.jf > AF100-12_5_reads.histo
jellyfish histo -t 10 SF2S9_reads.jf > SF2S9_reads.histo
jellyfish histo -t 10 AF100-12_7G_reads.jf > AF100-12_7G_reads.histo
