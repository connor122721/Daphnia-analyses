#!/bin/bash

#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH -t 8:00:00
#SBATCH --mem=64000
#SBATCH -p standard
#SBATCH -A berglandlab

chunk=$1

echo starting chunk ${chunk}

cd /scratch/kbb7sh/Daphnia/SingleMoms2018

echo mapping

module load gparallel/20170822
module load gcc/7.1.0
module load trimmomatic/0.36
module load bwa
module load samtools
module load picard

#Run program

parameterFile="/scratch/kbb7sh/Daphnia/SingleMoms2018/inputDoerthe"

varA=$(sed "${SLURM_ARRAY_TASK_ID}!d" $parameterFile | cut -f 1 )

bash trialmapping.sh $varA

