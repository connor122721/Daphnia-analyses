#!/usr/bin/env bash
#
#SBATCH -J rehh # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # 30 min
#SBATCH --mem 10G
#SBATCH -o /project/berglandlab/connor/rehh/rehh.out # Standard output
#SBATCH -e /project/berglandlab/connor/rehh/rehh.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load intel/18.0 intelmpi/18.0 R/3.6.3

# Working & temp directory
wd="/project/berglandlab/connor/rehh"

# Extract info
chrom=$( cat ${wd}/paramList | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# Run rehh
Rscript ${wd}/rehh.R ${chrom}
