#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/run_BLAST.log

#to remove all symlinks from the folder (and all subfolders)
#find . -type l -exec unlink {} \;

#to make symlinks to all files
cat genome_list_fullpath.txt | while read line 
do
   ln -s $line
done

mkdir genomes
mv *.fa genomes/
