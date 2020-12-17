#!/usr/bin/env bash
#
#SBATCH -J smc.joint # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-01:00 # 1 hours
#SBATCH --mem 5G
#SBATCH -o /scratch/csm6hg/smcpp/smcpp.out # Standard output
#SBATCH -e /scratch/csm6hg/smcpp/smcpp.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load SMC++ & Modules
module load gnuplot/5.2.2

# Working & temp directory
wd="/scratch/csm6hg/smcpp"
outdir="/scratch/csm6hg/smcpp/split"

# Extract input
samp1=$( cat ${wd}/samples_joint | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
samp2=$( cat ${wd}/samples_joint | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
chr=$( cat ${wd}/samples_joint | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
pop1=$( cat ${wd}/samples_joint | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
pop2=$( cat ${wd}/samples_joint | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f6 )

# Make conditional directories
[ ! -d ${outdir}/ ] && mkdir ${outdir}/
[ ! -d ${outdir}/${pop1}.${pop2}/ ] && mkdir ${outdir}/${pop1}.${pop2}/
[ ! -d ${outdir}/${pop2}.${pop1}/ ] && mkdir ${outdir}/${pop2}.${pop1}/

# Joint frequency spectrum SMC
./smc++ vcf2smc --cores 10 ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/final.sort.bed.gz ${outdir}/${pop1}.${pop2}/${pop1}.${pop2}.${chr}.pop12.smc.gz ${chr} ${pop1}:${samp1} ${pop2}:${samp2}
./smc++ vcf2smc --cores 10 ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/final.sort.bed.gz ${outdir}/${pop2}.${pop1}/${pop2}.${pop1}.${chr}.pop21.smc.gz ${chr} ${pop2}:${samp2} ${pop1}:${samp1}

# Finish
echo "Finish"