#!/usr/bin/env bash
#
#SBATCH -J smc++ # A single job name for the array
#SBATCH --ntasks-per-node=2 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 10G
#SBATCH -o /project/berglandlab/connor/smcpp/smcpp.out # Standard output
#SBATCH -e /project/berglandlab/connor/smcpp/smcpp.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

SLURM_ARRAY_TASK_ID=1

# Load SMC++ & Modules
source /home/csm6hg/smcpp/bin/activate
module load intel/18.0 intelmpi/18.0 R/3.6.3
module load gnuplot/5.2.2

# Working & temp directory
wd="/project/berglandlab/connor/smcpp"
outdir="/project/berglandlab/connor/smcpp/outdir"

# Extract input
samp=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
chr=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# VCF -> SMC
smc++ vcf2smc --cores 2 ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/final.sort.bed.gz ${outdir}/${samp}.${chr}.smc.gz ${chr} Pop1:${samp}

# Estimate demographic history
smc++ estimate --cores 2 5.69e-9 ${tmpdir}/${samp}.${chr}.smc.gz -o ${outdir}/${samp}/ --em-iterations 25

# Plot
smc++ plot --cores 2 -g 0.2 -c ${outdir}/${samp}/${samp}.${chr}.pdf ${outdir}/${samp}/*.json

# Remove temporary directory
rm -rf ${tmpdir}

# Finish
echo "Finish"
