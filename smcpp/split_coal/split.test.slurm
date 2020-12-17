#!/usr/bin/env bash
#
#SBATCH -J smc++ # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/csm6hg/smcpp/smcpp.out # Standard output
#SBATCH -e /scratch/csm6hg/smcpp/smcpp.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load SMC++ & Modules
module load gnuplot/5.2.2

# Working & temp directory
wd="/scratch/csm6hg/smcpp"
smc="/project/berglandlab/connor/smcpp/outdir"
outdir="/scratch/csm6hg/smcpp/split"

# Extract input
samp1=$( cat ${wd}/joint_split | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
samp2=$( cat ${wd}/joint_split | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
pop1=$( cat ${wd}/joint_split | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
pop2=$( cat ${wd}/joint_split | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )

# Estimate demographic history
./smc++ split -o ${outdir}/split/ --base ${pop1}_${pop2}_split --em-iterations 30 --timepoints 1 1e8 ${smc}/${samp1}/model.final.json ${smc}/${samp2}/model.final.json ${smc}/smc/${samp1}.*.smc.gz ${smc}/smc/${samp2}.*.smc.gz ${outdir}/${pop1}_${pop2}/${pop1}.${pop2}*smc.gz ${outdir}/${pop2}_${pop1}/${pop2}.${pop1}*smc.gz

# Plot
./smc++ plot --cores 10 -g 0.2 -c ${outdir}/split/${pop1}_${pop2}.pdf ${outdir}/split/${pop1}*${pop2}*split.final.json

# Finish
echo "Finish"
