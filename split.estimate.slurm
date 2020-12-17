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

SLURM_ARRAY_TASK_ID=1

# Load SMC++ & Modules
source /home/csm6hg/smcpp/bin/activate
module load gnuplot/5.2.2

# Working & temp directory
wd="/scratch/csm6hg/smcpp"
outdir="/scratch/csm6hg/smcpp/split"

# Extract input
samp1=$( cat ${wd}/remaining_splits | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
samp2=$( cat ${wd}/remaining_splits | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )
pop1=$( cat ${wd}/remaining_splits | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )
pop2=$( cat ${wd}/remaining_splits | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f5 )
list="Scaffold_1863_HRSCAF_2081 Scaffold_1931_HRSCAF_2197 Scaffold_2158_HRSCAF_2565 Scaffold_2217_HRSCAF_2652 Scaffold_2373_HRSCAF_2879 Scaffold_6786_HRSCAF_7541 Scaffold_7757_HRSCAF_8726 Scaffold_9197_HRSCAF_10753 Scaffold_9198_HRSCAF_10754 Scaffold_9199_HRSCAF_10755 Scaffold_9200_HRSCAF_10757 Scaffold_9201_HRSCAF_10758"

# Make conditional directories
[ ! -d ${outdir}/ ] && mkdir ${outdir}/
[ ! -d ${outdir}/${pop1}_${pop2}/ ] && mkdir ${outdir}/${pop1}_${pop2}/
[ ! -d ${outdir}/${pop2}_${pop1}/ ] && mkdir ${outdir}/${pop2}_${pop1}/

# Forloop through chromosome
for i in $list; do
smc++ vcf2smc --cores 10 ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/final.sort.bed.gz ${outdir}/${pop1}_${pop2}/${pop1}.${pop2}.${i}.pop12.smc.gz ${i} ${samp1}:${samp1} ${samp2}:${samp2}
smc++ vcf2smc --cores 10 ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/final.sort.bed.gz ${outdir}/${pop2}_${pop1}/${pop2}.${pop1}.${i}.pop21.smc.gz ${i} ${samp2}:${samp2} ${samp1}:${samp1}
done

# WHole genome of each individual
smc++ estimate --base ${pop1}_split --em-iterations 30 --timepoints 1 1e8 --cores 10 5.69e-09 ${outdir}/${pop1}_${pop2}/*.smc.gz
smc++ estimate--base ${pop2}_split --em-iterations 30 --timepoints 1 1e8 --cores 10 5.69e-09 ${outdir}/${pop2}_${pop1}/*.smc.gz

# Estimate demographic history
smc++ split -o ${outdir}/split/ --base ${pop1}_${pop2}_split --em-iterations 30 --timepoints 1 1e8 --cores 10 ${smc}/${samp1}/model.final.json ${smc}/${samp2}/model.final.json ${smc}/smc/${samp1}*.smc.gz ${smc}/smc/${samp2}*.smc.gz ${outdir}/${pop1}_${pop2}/${pop1}*${pop2}*.smc.gz ${outdir}/${pop2}_${pop1}/${pop2}*${pop1}*.smc.gz

# Plot
smc++ plot --cores 10 -g 0.2 -c ${outdir}/split/${pop1}_${pop2}.pdf ${outdir}/split/${pop1}*${pop2}*split.final.json

# Finish
echo "Finish"
