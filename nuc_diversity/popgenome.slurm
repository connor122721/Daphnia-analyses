#!/usr/bin/env bash
#
#SBATCH -J slim-constant # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:30 # 30 min
#SBATCH --mem 15G
#SBATCH -o /project/berglandlab/connor/PopGenome/scripts/SC.pop.out # Standard output
#SBATCH -e /project/berglandlab/connor/PopGenome/scripts/SC.pop.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load intel/18.0 intelmpi/18.0 R/3.6.3

# Working & temp directory
wd="/project/berglandlab/connor/PopGenome/scripts"

# Extract constants from parameter file
SC=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
Chrom=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )
# Rep=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f8 )
# Samp=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
UpMAF=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
LoMAF=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
Out=$( cat ${wd}/SC_model_paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# RAM disk
[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/ 
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID} 
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

# Run PopGenome
Rscript ${wd}/diversity.SC.updated.R ${SC} ${Chrom} ${UpMAF} ${LoMAF} ${Out} ${tmpdir}

# Remove temp directory
rm -rf tmpdir
