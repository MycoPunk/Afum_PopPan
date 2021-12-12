#!/usr/bin/bash -l
#SBATCH -p short -C xeon -N 1 -n 96 --mem 24gb --out logs/orthofinder.%A.log

mkdir -p logs
module load orthofinder/2.5.2
opt="-C xeon" # could change to "-C xeon" and will run on the xeon nodes; # could change this to empty and will run on any node
JOBS=orthofinder_steps.diamond.sh
LOG=orthofinder_steps.diamond.log
CHUNK=50
export TMPDIR=/scratch
if [ ! -f $LOG ]; then
	orthofinder -op -t 96 -a 96 -f input -S diamond_ultra_sens -o OrthoFinder_diamond > $LOG
fi
grep ^diamond $LOG | grep -v 'commands that must be run' | perl -p -e 's/-p 1/-p 8/g'> $JOBS

t=$(wc -l $JOBS | awk '{print $1}')
MAX=$(expr $t / $CHUNK)
echo "t is $t MAX is $MAX"
for n in $(seq $MAX)
do
	START=$(perl -e "printf('%d',1 + $CHUNK * ($n - 1))")
	END=$(perl -e "printf('%d',$CHUNK* $n)")
#	echo "$START,$END for $n"
	run=$(sed -n ${START},${END}p $JOBS)
#		echo "sbatch -p short -N 1 -n 1 --mem 4gb --wrap \"module load orthofinder/2.5.2; $run\""
	sbatch $opt --out logs/diamond.$n.log -J Dmd$n -p short -N 1 -n 8 --mem 4gb --wrap "module load orthofinder/2.5.2; $run"
done
