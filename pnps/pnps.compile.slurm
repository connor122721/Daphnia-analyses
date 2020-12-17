#!/usr/bin/env bash
#
#SBATCH -J pnps # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:10 # 10 min
#SBATCH --mem 5G
#SBATCH -o /project/berglandlab/connor/pnps/scripts/pnps.out # Standard output
#SBATCH -e /project/berglandlab/connor/pnps/scripts/pnps.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load intel/18.0 intelmpi/18.0 R/3.6.3

# Working & temp directory
wd="/project/berglandlab/connor/pnps/scripts"
outdir="/project/berglandlab/connor/pnps/data/output"

# Run PopGenome
Rscript ${wd}/compile.write.pnps.R

# Remove files
rm -rf ${outdir}

# Finish
echo "Finish"