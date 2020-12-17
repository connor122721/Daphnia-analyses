#!/usr/bin/env bash
#
#SBATCH -J smc++ # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 7G
#SBATCH -o /scratch/csm6hg/smcpp/smcpp.out # Standard output
#SBATCH -e /scratch/csm6hg/smcpp/smcpp.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

SLURM_ARRAY_TASK_ID=1

# Load SMC++ & Modules
source /home/csm6hg/smcpp/bin/activate
module load intel/18.0 intelmpi/18.0 R/3.6.3
module load gnuplot/5.2.2

# Working & temp directory
wd="/scratch/csm6hg/smcpp"
smc="/project/berglandlab/connor/smcpp/outdir"
outdir="/scratch/csm6hg/smcpp/split"

# Extract input
samp1=$( cat ${wd}/joint_split | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
samp2=$( cat ${wd}/joint_split | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# Estimate demographic history
smc++ split 
-o ${outdir}/split/ 
--base split.${samp1}.${samp2} 
--cores 10 
--em-iterations 30 
--timepoints 1 1e8 
${smc}/${samp1}/model.final.json 
${smc}/${samp2}/model.final.json 
${smc}/smc/${samp1}*.smc.gz 
${smc}/smc/${samp2}*.smc.gz 
${outdir}/${samp1}.${samp2}.*.smc.gz 
${outdir}/${samp2}.${samp1}.*.smc.gz

# Plot
smc++ plot --cores 5 -g 0.2 -c ${outdir}/${samp}.pdf ${outdir}/split/split.${samp1}.${samp2}*

# Finish
echo "Finish"