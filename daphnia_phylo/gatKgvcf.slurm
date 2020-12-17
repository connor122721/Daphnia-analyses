#!/usr/bin/env bash
#
#SBATCH -J gatk # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 120G
#SBATCH -o /scratch/csm6hg/daphnia_phylo/gatk.out # Standard output
#SBATCH -e /scratch/csm6hg/daphnia_phylo/gatk.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load modules
module load gatk/4.0.0.0

# Parameters
parameterFile="/scratch/csm6hg/daphnia_phylo/Daphnia.SRR.csv"
wd="/scratch/csm6hg"

# Extract sample
samp=$( cat ${parameterFile} | grep "^$SLURM_ARRAY_TASK_ID," | cut -d',' -f4 )

# Add read groups
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
-I ${wd}/bam/${samp}_finalmap.bam \
-O ${wd}/bam/.bam \
-LB "library" \
-PL "ILLumina" \
-PU "platunit" \
-SM ${samp} 

# Run program
gatk --java-options "-Xmx4g" HaplotypeCaller --native-pair-hmm-threads 10 \
-R ${wd}/daphnia_phylo/totalHiCwithallbestgapclosed.fa \
-I ${wd}/bam/${samp}_finalmap.bam \
-O ${wd}/daphnia_phylo/gvcf/${samp}.g.vcf \
-ERC GVCF 