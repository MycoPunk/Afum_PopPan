#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 --mem 4gb
cd /bigdata/stajichlab/shared/projects/Afumigatus_pangenome/Pangenome_analysis/PIRATE/260_strains

#this grabs the family names
gene_fams=$(grep '^>' representative_sequences.faa | sed -e 's/;.*//' | sed 's/>//1')

#this grabs the mRNA name
mRNA=$(grep 'locus_tag=' representative_sequences.faa | sed 's/^.*locus_tag=//'| sed -e 's/;.*//')

#this slaps them together
paste <(printf %s "$gene_fams") <(printf %s "$mRNA") >rep_gene_fams_w_mRNA.txt


#get the positions of representative gene families relative to each assembly
names=$(awk '{print $2}' rep_gene_fams_w_mRNA.txt)
pos=$(co-ords/*)

while read LINE; do
	grep -f < "$LINE" "${pos}"
done <"${names}"


cat rep_gene_fams_w_mRNA.txt | while read ORTHOGENEID MRNAID
do
 	STRAIN=$(echo $MRNAID | perl -p -e 's/_mRNA_\S+$//')
	echo -e -n "$ORTHOGENEID\t$MRNAID\t"; 
	grep $MRNAID co-ords/${STRAIN}.co-ords.tab | cut -f2,3,4,8 
done > rep_gene_fams_w_mRNA.with_coords.tsv


####get denovo  AF293 gene fams for validation
#grep -H "AF293" * >AF293_gene_fams.txt

#sed -i 's/.nucleotide.fasta:>/ /g' AF293_gene_fams.txt
