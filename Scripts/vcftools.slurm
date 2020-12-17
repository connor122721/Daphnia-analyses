#!/usr/bin/env bash
#
#
#SBATCH -J vcftools-job1 # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-03:00 # Running time of 15 minutes
#SBATCH --mem 5G # Memory request of 2 GB
#SBATCH -o /home/project/berglandlab/connor/PopGenome/vcf/vcftools.out # Standard output
#SBATCH -e /home/project/berglandlab/connor/PopGenome/vcf/vcftools.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Modules
module load vcftools/0.1.15
module load htslib/1.9

# Removes all sites > 0.25 MAF from filtered vcf.gz
vcftools --gzvcf MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.vcf.gz --recode --out maf_filter_0.25.ann.vcf.gz --max-maf 0.25 --chr goodChrom.txt --snps finalsetsnpset01pulex_table_wpulicaria_20200401

# gzip and tabix VCF
bgzip maf_filter_goodpos.ann.vcf.gz.recode.vcf > maf_filter_goodpos.recode.vcf.gz
tabix -p vcf maf_filter_0.25.out.recode.vcf.gz

# Theta pi
vcftools --gzvcf maf_filter_0.25.out.recode.vcf.gz --chr goodChrom.txt --window-pi 10000 --out maf_filter_0.25_pi

# Tajima’s D
vcftools --gzvcf maf_filter_0.25.out.recode.vcf.gz --chr goodChrom.txt --TajimaD 10000 --out maf_filter_0.25_tajimasD


vcftools --gzvcf maf_filter_0.25.recode.vcf.gz --recode --out maf_filter_goodpos.ann.vcf.gz --positions goodPOS.txt
