#!/usr/bin/env bash
#
#SBATCH -J smc++ # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 5G
#SBATCH -o /project/berglandlab/connor/smcpp/output/smcpp.out # Standard output
#SBATCH -e /project/berglandlab/connor/smcpp/output/smcpp.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load SMC++ & Modules
# source /home/csm6hg/smcpp/bin/activate
module load intel/18.0 intelmpi/18.0 R/3.6.3
module load gnuplot/5.2.2

SLURM_ARRAY_TASK_ID=3
# Working & temp directory
wd="/project/berglandlab/connor/smcpp"
outdir="/project/berglandlab/connor/smcpp/output"

# RAM disk
[ ! -d /dev/shm/$USER/ ] && mkdir /dev/shm/$USER/
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}
[ ! -d /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} ] && mkdir /dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

tmpdir=/dev/shm/$USER/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

# Extract input
samp=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
chr=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# VCF -> SMC
# smc++ vcf2smc ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/final.bed.gz ${wd}/${samp}.${chr}.smc.gz ${chr} Pop1:${samp}
./smc++ vcf2smc ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/NsandDepthandChrEnd.sorted.500merged.bed.gz ${wd}/A.${chr}.smc.gz ${chr} A:April_2017_D8_103,April_2017_D8_119
./smc++ vcf2smc ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/NsandDepthandChrEnd.sorted.500merged.bed.gz ${wd}/D.${chr}.smc.gz ${chr} D:Fall_2016_D10_63,Fall_2016_D10_67

# Estimate demographic history
# smc++ estimate 5.69e-9 ${wd}/${samp}.${chr}.smc.gz -o ${outdir}
./smc++ estimate 5.69e-9 ${wd}/A.${chr}.smc.gz -o ${outdir}/A
./smc++ estimate 5.69e-9 ${wd}/D.${chr}.smc.gz -o ${outdir}/D

# Plot and output data
smc++ plot ${wd}/${samp}.${chr}.smcpp.pdf ${outdir}/*.json -c -g 5
./smc++ plot ${wd}/A.smcpp.test.pdf ${outdir}/A/*.json -c -g 0.2 -x 5e3 4e4 -y 3e4 2e5
./smc++ plot ${wd}/D.smcpp.test.pdf ${outdir}/D/*.json -c -g 0.2 -x 5e3 4e4 -y 3e4 2e5

# Joint freqeuncy spectrum
./smc++ vcf2smc ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/NsandDepthandChrEnd.sorted.500merged.bed.gz ${wd}/A.joint.${chr}.smc.gz ${chr} A:April_2017_D8_103,April_2017_D8_119 D:Fall_2016_D10_63,Fall_2016_D10_67
./smc++ vcf2smc ${wd}/MapJune2020_ann.vcf.gz -m ${wd}/NsandDepthandChrEnd.sorted.500merged.bed.gz ${wd}/D.joint.${chr}.smc.gz ${chr} D:Fall_2016_D10_63,Fall_2016_D10_67 A:April_2017_D8_103,April_2017_D8_119

# Split
./smc++ split -o ${outdir}/ ${outdir}/A/model.final.json ${outdir}/D/model.final.json ${wd}/*.joint.*.smc.gz --em-iterations 20
./smc++ plot joint.pdf ${outdir}/model.final.json -c -g 0.2 -x 5e3 4e4 -y 3e4 2e5

# Remove tmpdir
rm -rf tmpdir

# Finish
echo "Finish"
