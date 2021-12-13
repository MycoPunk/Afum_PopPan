#!/usr/bin/bash
#SBATCH -p short --mem 32gb -N 1 -n 32 --out logs/asm_map.log

module load minimap2
module unload perl
module load parallel
URL1=https://fungidb.org/common/downloads/release-46/AfumigatusA1163/fasta/data/FungiDB-46_AfumigatusA1163_Genome.fasta
URL2=https://fungidb.org/common/downloads/release-46/AfumigatusAf293/fasta/data/FungiDB-46_AfumigatusAf293_Genome.fasta
REFFOLDER=ref_genome
mkdir -p $REFFOLDER
pushd $REFFOLDER
for n in $URL1 $URL2
do
	f=$(basename $n)
	if [ ! -s $f ]; then
		curl -O $n
		minimap2 -x asm10 -d $n.mmi $n
	fi
done
popd
#use Af293 for now
REFGENOME=$REFFOLDER/$(basename $URL2)
REFNAME=Af293
OUTDIR=asm_mapping
CLUSTER_MAP=cluster_mapping

mkdir -p $OUTDIR $CLUSTER_MAP

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes
ANNOTATEDIR=$(realpath annotate)
SAMPFILE=samples.csv

parallel -j $CPU minimap2 -t 2 --cs -cx asm10 $REFGENOME $INDIR/{}.masked.fasta \> $OUTDIR/{}.$REFNAME.paf ::: $(cut -f1 -d, $SAMPFILE)
parallel -j $CPU ln -s $ANNOTATEDIR/{}/annotate_misc/antismash/clusters.bed $CLUSTER_MAP/{}.clusters.bed ::: $(cut -f1 -d, $SAMPFILE)

parallel -j $CPU paftools.js liftover -l 5000 $OUTDIR/{}.$REFNAME.paf $CLUSTER_MAP/{}.clusters.bed \> $CLUSTER_MAP/{}.clusters_mapped_Af293.bed ::: $(cut -f1 -d, $SAMPFILE)

module load bedtools
cat cluster_mapping/*.clusters_mapped_Af293.bed | bedtools sort -i - | bedtools merge -c 1 -o count -i - -d 5000  > cluster_frequency.tab

cat cluster_mapping/*.clusters_mapped_Af293.bed | bedtools sort -i - | bedtools cluster -i - -d 5000 > cluster_names.tab
