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
module load samtools/1.10
module load bcftools/1.9
module load htslib/1.9
module load psmc/0.6.5
module load gnuplot/5.2.2


# Working & temp directory
wd="/project/berglandlab/connor/psmc/output"
outdir="/project/berglandlab/connor/psmc/output/vcf"
bam="/project/berglandlab/Karen/MappingDec2019/bams/PulexBams"

SLURM_ARRAY_TASK_ID=2100
# Extract input
samp=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f2 )
chr=$( cat ${wd}/sample_list | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f3 )

# BAM to VCF.gz
bcftools mpileup -f ${wd}/totalHiCwithallbestgapclosed.fa ${bam}/${samp}_finalmap_mdup.bam -t ${chr} -T ${wd}/final.bed.gz -t ${chr} | bcftools call -cv -Oz -o ${outdir}/${samp}.${chr}.vcf.gz
bcftools index ${outdir}/${samp}.${chr}.vcf.gz -t
bcftools consensus -f ${wd}/totalHiCwithallbestgapclosed.fa ${outdir}/${samp}.${chr}.vcf.gz | vcfutils.pl vcf2fq > test.fq

# VCF to FQ
vcfutils.pl vcf2fq ${outdir}/${samp}.${chr}.vcf.gz > test.fq

# Combine FQ
# cat ${outdir}/${samp}.${chr}.fq > ${outdir}/${samp}.fq

# FQ -> PSMCFA
#/project/berglandlab/connor/psmc/utils/fq2psmcfa -q 20 test.fq.gz > ${outdir}/${samp}.psmcfa

# PSMC
#psmc -N25 -t15 -r1.75 -p "4+25*2+4+6" -o ${outdir}/${samp}.psmc ${outdir}/${samp}.psmcfa

# Create pdf psmc plot
#../utils/psmc_plot.pl -u 5.69e-09 -g 0.05 ${samp} ${wd}/${samp}.psmc 

# Finish
echo "Finish"
