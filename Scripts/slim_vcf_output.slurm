#!/usr/bin/env bash
#
#
#SBATCH -J slim-vcf-output # A single job name for the array
#SBATCH --ntasks-per-node=1 # one core
#SBATCH -N 1 # on one node
#SBATCH --cpus-per-task=1 ### standard has 28 or 40 $SLURM_CPUS_PER_TASK
#SBATCH -t 0-04:00 # Running time of 4hr
#SBATCH --mem 10G # Memory request of 10 GB
#SBATCH -o /project/berglandlab/connor/slim/scripts/slim.vcf.output.out # Standard output
#SBATCH -e /project/berglandlab/connor/slim/scripts/slim.vcf.output.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

### run command
module load samtools/1.10
module load htslib/1.9
module load intel/18.0 intelmpi/18.0 R/3.6.3

# Change directory
DIR=$(constant6)

# Make new directory
mkdir /project/berglandlab/connor/slim/slim_output/${DIR}

# Move vcf files
mv /project/berglandlab/connor/slim/slim_output/*.vcf /project/berglandlab/connor/slim/slim_output/${DIR}

# gzip and tabix vcf files
find /project/berglandlab/connor/slim/slim_output/constant75-100 -maxdepth 1 -type f -exec bgzip {} \;
find /project/berglandlab/connor/slim/slim_output/constant75-100 -maxdepth 1 -type f -exec tabix {} \;

# Diversity script [data.path] [output]
Rscript /project/berglandlab/connor/slim/scripts/slim_vcf_output.R /project/berglandlab/connor/slim/slim_output/${DIR}/ ${DIR}.output.csv

# Population script [data.path] [output]
Rscript /project/berglandlab/connor/slim/scripts/compile.write.popsize.R /project/berglandlab/connor/slim/slim_output/${DIR}/ ${DIR}.pop.csv

# Move population files
mv /project/berglandlab/connor/slim/slim_output/*.vcf /project/berglandlab/connor/slim/slim_output/${DIR}

