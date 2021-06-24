#!/usr/bin/bash
#SBATCH --mem 10G --ntasks 20 --nodes 1  --time 2:00:00 

cd /bigdata/stajichlab/shared/projects/Afumigatus_pangenome/annotate

#make a copy and name each file with the folder two above (genome name)
for subdir in *; do 
	cp "$subdir"/annotate_misc/annotations.dbCAN.txt "$subdir"/annotate_misc/"$subdir".CAZY.txt; 
done;

#move the copied files into a new dir
for subdir in *; do 
	cp "$subdir"/annotate_misc/*CAZY.txt ~/bigdata/Afum_pan_genome/CAZY/; 
done;


#ok, so there are 315 CAZY sets, and we only want the ones in our pan genome. 
#import strain list used in pirate. 
cd ~/bigdata/Afum_pan_genome/CAZY
cp /bigdata/stajichlab/shared/projects/Afumigatus_pangenome/Pangenome_analysis/PIRATE/261_strains/genome_list.txt .

#clean up the file names so they match 
sed -i 's/Afum_//g' genome_list.txt
sed -i 's/$/.CAZY.txt/' genome_list.txt

#clean up the files
#replace all "-" with "_"

find . -type f -name '*.txt' | while read FILE ; do
    newfile="$(echo ${FILE} |sed -e 's/\-/_/g')" ;
    mv "${FILE}" "${newfile}" ;
done

#replace all "." before .CAZY with "_"
for f in *.*.*.txt; do 
	i="${f%.CAZY.txt}"; mv -i -- "$f" "${i//./_}.CAZY.txt"; 
done


#copy the matches into a new folder
mkdir CAZY_261

while read file; do 
	cp "$file" CAZY_261/; 
done < genome_list.txt

#check that you're not missing any
ls | wc -l
#261 #nice. 
