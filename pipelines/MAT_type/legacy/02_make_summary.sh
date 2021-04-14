#!/usr/bin/bash
#SBATCH --mem=20G -p batch --nodes 1 --ntasks 2 --out logs/make_summary.log

#this script writes a summary file from the blast results for mating type genes

cd results
MAT1=""
NAME=""
for i in *MAT1-1.out; do
name=$(basename "$i" .scaffolds.fa.rename.fa.MAT1-1.out)
val=$(awk 'FNR == 6 {print $12}' $i)
  NAME+="${name}\n"	
  MAT1+="${val}\n"
done

MAT2=""
for i in *MAT1-2.out; do
val=$(awk 'FNR == 6 {print $12}' $i)
  MAT2+="${val}\n"
done

echo -e "$(paste -d '\t' <(echo -e "$NAME") <(echo -e "$MAT1") <(echo -e "$MAT2"))" >results_summary1.tab

#compare 
awk '{
if ($2 > $3)
	print $1, $2, $3, "MAT1";
else
	print $1, $2, $3, "MAT2";
}' results_summary1.tab > results_summary2.tab

rm results_summary1.tab

#add headers
sed -i 1i"strain,MAT1-1_bitscore,MAT1-2_bitscore, MAT_type" results_summary2.tab


####
#clean up the summary file
awk '{gsub("Aspergillus_fumigatus_", "");
        print;  
}' results_summary2.tab > results_summary3.tab

rm results_summary2.tab

#remove that trailing line
tac results_summary3.tab | sed "1,1d" | tac > MATresults_summary.tab

rm results_summary3.tab
