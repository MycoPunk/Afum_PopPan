#!/bin/bash
#SBATCH --nodes 1 --ntasks 24 --mem 192gb -J afumAAFTF --out logs/AAFTF_full.%a.%A.log -p batch --time 72:00:00

hostname
MEM=192
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

module load AAFTF

OUTDIR=input
SAMPLEFILE=samples.dat
ASM=genomes
WORKDIR=working_AAFTF
mkdir -p $ASM
mkdir -p $WORKDIR
if [ -z $CPU ]; then
    CPU=1
fi

sed -n ${N}p $SAMPLEFILE | while read BASE STRAIN PHYLUM
do

    ASMFILE=$ASM/${STRAIN}.spades.fasta

    VECCLEAN=$ASM/${STRAIN}.vecscreen.fasta
    PURGE=$ASM/${STRAIN}.sourpurge.fasta
    CLEANDUP=$ASM/${STRAIN}.rmdup.fasta
    PILON=$ASM/${STRAIN}.pilon.fasta
    SORTED=$ASM/${STRAIN}.sorted.fasta
    STATS=$ASM/${STRAIN}.sorted.stats.txt
    LEFTTRIM=$WORKDIR/${BASE}_1P.fastq.gz
    RIGHTTRIM=$WORKDIR/${BASE}_2P.fastq.gz

    LEFT=$WORKDIR/${BASE}_filtered_1.fastq.gz
    RIGHT=$WORKDIR/${BASE}_filtered_2.fastq.gz

    echo "$BASE $STRAIN"
    if [ ! -f $ASMFILE ]; then    
	if [ ! -f $LEFT ]; then
	    echo "$OUTDIR/${BASE}_R1.fq.gz $OUTDIR/${BASE}_R2.fq.gz"
	    if [ ! -f $LEFTTRIM ]; then
		if [[ $BASE = "ERR062316" || $BASE = "ERR062317" ]]; then
		    AAFTF trim --method bbduk --memory $MEM --left $OUTDIR/${BASE}_R1.fq.gz --right $OUTDIR/${BASE}_R2.fq.gz -c $CPU -o $WORKDIR/${BASE} --minlength 35
		else
		    AAFTF trim --method bbduk --memory $MEM --left $OUTDIR/${BASE}_R1.fq.gz --right $OUTDIR/${BASE}_R2.fq.gz -c $CPU -o $WORKDIR/${BASE}
		fi
	    fi
	    echo "$LEFTTRIM $RIGHTTRIM"
	    AAFTF filter -c $CPU --memory $MEM -o $WORKDIR/${BASE} --left $LEFTTRIM --right $RIGHTTRIM --aligner bbduk 
	    echo "$LEFT $RIGHT"
	    if [ -f $LEFT ]; then
		unlink $LEFTTRIM
		unlink $RIGHTTRIM
	    fi
	fi
	AAFTF assemble -c $CPU --left $LEFT --right $RIGHT  --memory $MEM \
	    -o $ASMFILE -w $WORKDIR/spades_${STRAIN}
	
	if [ -s $ASMFILE ]; then
	    rm -rf $WORKDIR/spades_${STRAIN}/K?? $WORKDIR/spades_${STRAIN}/tmp
	fi
	
	if [ ! -f $ASMFILE ]; then
	    echo "SPADES must have failed, exiting"
	    exit
	fi
    fi
    
    if [ ! -f $VECCLEAN ]; then
	AAFTF vecscreen -i $ASMFILE -c $CPU -o $VECCLEAN 
    fi

    if [ ! -f $PURGE ]; then
	AAFTF sourpurge -i $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM --left $LEFT  --right $RIGHT
    fi
    
    if [ ! -f $CLEANDUP ]; then
	AAFTF rmdup -i $PURGE -o $CLEANDUP -c $CPU -m 500
    fi
    
    if [ ! -f $PILON ]; then
	AAFTF pilon -i $CLEANDUP -o $PILON -c $CPU --left $LEFT  --right $RIGHT
    fi
    
    if [ ! -f $PILON ]; then
	echo "Error running Pilon, did not create file. Exiting"
	exit
    fi
    
    if [ ! -f $SORTED ]; then
	AAFTF sort -i $PILON -o $SORTED
    fi
    
    if [ ! -f $STATS ]; then
	AAFTF assess -i $SORTED -r $STATS
    fi
done
