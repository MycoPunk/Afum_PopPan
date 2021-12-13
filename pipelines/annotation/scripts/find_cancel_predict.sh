echo -n "sbatch -p batch -a "
for n in $(ls logs/predict.*.log); do ct=$(grep -c "CANCEL" $n); if [ "$ct" -gt 0 ]; then  echo $n | perl -p -e 's/\S+\.(\d+)\.log/$1/'; fi; done | perl -p -e 's/\n/,/' && echo " pipeline/02_predict.sh"
