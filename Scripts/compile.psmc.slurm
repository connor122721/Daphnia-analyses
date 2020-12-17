#!/usr/bin/env bash
#
#SBATCH -J psmc.comp # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-00:10 # 10 min
#SBATCH --mem 1G
#SBATCH -o /project/berglandlab/connor/psmc/output/compile.out # Standard output
#SBATCH -e /project/berglandlab/connor/psmc/output/compile.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load intel/18.0 intelmpi/18.0 R/3.6.3

# Working directory
wd="/project/berglandlab/connor/psmc/output"
file="/project/berglandlab/connor/psmc/output/vcf"

# Compile psmc results: files; output; mutation; segments length; gen/year
Rscript ${wd}/compile.write.psmc.R  ${file} psmc.out.csv 5.69e-09 100 0.2

# Finish
echo "Finish"
