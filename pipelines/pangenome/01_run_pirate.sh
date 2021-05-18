#!/usr/bin/bash
#SBATCH -p batch -N 1 -n 32 --mem 128gb  --out logs/optimized_261_strains.log -J PIRATE

module unload miniconda3
module unload miniconda2
module unload python
module unload perl
module unload anaconda2
module load anaconda3
module load cd-hit/4.8.1
module load mcl
module load ncbi-blast/2.9.0+
module load diamond
export TEMPDIR=/scratch
export TMPDIR=/scratch
CPUS=$SLURM_CPUS_ON_NODE
if [ -z $CPUS ]; then
 CPUS=1
fi

#source activate pirate
source activate /bigdata/stajichlab/jstajich/.conda/envs/pirate

#~/projects/PIRATE/bin/PIRATE -f "mRNA" -i gff \
/bigdata/stajichlab/jstajich/projects/PIRATE/bin/PIRATE -f "mRNA" -i gff \
    -t $CPUS --rplots -o 261_strains --nucl \
    -s "85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100" -k "--cd-low 100 -e 1E-9 --hsp-prop 0.5" -a -r
