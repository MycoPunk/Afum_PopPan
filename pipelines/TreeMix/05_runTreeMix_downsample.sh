#!/bin/bash
#SBATCH --job-name=run_TreeMix_downwample
#SBATCH --output=logs/run_TreeMix_downsample.%a.log
#SBATCH --mail-user=lotuslofgren@gmail.com
#SBATCH --mail-type=FAIL,END
#SBATCH --time=4:00:00
#SBATCH --mem=15G
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=10
#SBATCH --partition=common

module load TreeMix/1.13


for m in {1..10}; do
    for i in {1..10}; do
        # Generate random seed
        s=$RANDOM
	treemix -i VCF_samp_${i}.recode_noinv.treemix.frq.gz -o A_fum_subpop.${i}.${m} -global -m ${m} -seed ${s}
    done
done


#pack them up
mkdir TreeMix_results_OG
mv A_fum_subpop* TreeMix_results_OG
#tar -zcvf treemix.tar.gz TreeMix_results
