#!/usr/bin/env bash
#
#SBATCH -J psmc.test # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 7G
#SBATCH -o /project/berglandlab/connor/psmc/output/psmc.test.out # Standard output
#SBATCH -e /project/berglandlab/connor/psmc/output/psmc.test.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load psmc/0.6.5
module load gnuplot/5.2.2
module load ghostscript/9.52

SLURM_ARRAY_TASK_ID=264

# Working & temp directory
wd="/project/berglandlab/connor/psmc/output"
outdir="/project/berglandlab/connor/psmc/output/vcf"
bam="/project/berglandlab/Karen/MappingDec2019/bams/PulexBams"
home="/home/csm6hg/samtools-0.1.19"

# Extract input
samp=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )

# BAM to VCF.gz
cd ${home}
./samtools mpileup -C50 -uf ${wd}/totalHiCwithallbestgapclosed.fa ${bam}/${samp}_finalmap_mdup.bam -l ${wd}/final.bed.gz | ./bcftools/bcftools view -c - \ | ./bcftools/vcfutils.pl vcf2fq -d 10 -D 100 -Q 30 | gzip > ${samp}.fq.gz

# FQ -> PSMCFA
fq2psmcfa ${outdir}/${samp}.fq.gz > ${outdir}/${samp}.psmcfa

# PSMC
psmc -N25 -t15 -r1.75 -p "4+25*2+4+6" -d -o ${outdir}/${samp}.psmc ${outdir}/${samp}.psmcfa

# Create pdf psmc plot
/project/berglandlab/connor/psmc/utils/psmc_plot.pl -p -u 5.69e-09 -g 0.2 ${outdir}/${samp} ${outdir}/${samp}.psmc 

# Finish
echo "Finish"
