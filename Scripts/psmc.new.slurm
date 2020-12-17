#!/usr/bin/env bash
#
#SBATCH -J msmc.test # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 5G
#SBATCH -o /project/berglandlab/connor/msmc-tools/output/msms.out # Standard output
#SBATCH -e /project/berglandlab/connor/msmc-tools/output/msmc.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gcc/9.2.0 mvapich2/2.3.3
module load python/3.7.7

# Working & temp directory
wd="/project/berglandlab/connor/msmc-tools/output"
outdir="/project/berglandlab/connor/msmc-tools/output/vcf"

# Extract input
samp=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
chr=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Generate mutlihet file
./generate_multihetsep.py --mask ${outdir}/final.bed.gz --chr ${chr} ${outdir}/${samp}.${chr}.vcf.gz

# Finish
echo "Finish"
