
for file in in $(cat poor_samples.csv); do if [ -f genomes/$file.scaffolds.fasta ]; then  echo "mv genomes/$file.scaffolds.fasta bad_genomes"; fi; done
