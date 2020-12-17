#!/usr/bin/env bash
#
#SBATCH -J vcftools # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-02:00 # 2 hours
#SBATCH --mem 10G
#SBATCH -o /project/berglandlab/connor/smcpp/vcf.out # Standard output
#SBATCH -e /project/berglandlab/connor/smcpp/vcf.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gcc/9.2.0
module load bedtools/2.29.2
module load htslib/1.10.2
module load vcftools/0.1.16

# Remove multiple enteries
vcftools --gzvcf MapJune2020_ann.vcf.gz --max-alleles 2 --recode --out MapJune2020_ann_intersect_mult
bgzip MapJune2020_ann_intersect_mult.vcf
tabix -p MapJune2020_ann_intersect_mult.vcf