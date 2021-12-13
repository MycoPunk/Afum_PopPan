
for n in $(ls *.sorted.fasta); do b=$(basename $n .sorted.fasta); re=$(grep "^$b" ../samples.csv);  if [ -z "$re" ]; then new=$(echo $b | perl -p -e 's/_/-/g'); echo "$b,$new"; fi; done > ../samples_new
