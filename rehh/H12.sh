#! /bin/bash


#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G
#SBATCH --time=2:00:00
#SBATCH --partition=standard
#SBATCH --account=bergland-erickson
#SBATCH --array=1-40

intervals=/scratch/pae3g/zaps/ref/z_indianus_16GNV01_v02.intervals

i=`sed -n ${SLURM_ARRAY_TASK_ID}p $intervals`

cd /scratch/pae3g/zaps/popgen

module load bcftools
#this pulls out genotypes and subsets samples to unrelated individuals, removing any invariant sites with --min-ac=1
bcftools   annotate -x INFO,^FORMAT/GT /scratch/pae3g/zaps/haplotype_calling/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.${i}.recode.vcf.gz -Ou |  bcftools view -S /scratch/pae3g/zaps/popgen/CM_HPO_unrelated_ids.txt --min-ac=1 -Ov -o /scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.${i}.GT_unrelatedCMHPO.vcf

module load gcc/7.1.0  mvapich2/2.3.1 python/2.7.16

python /scratch/pae3g/zaps/scripts/makeH12file.py /scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.${i}.GT_unrelatedCMHPO.vcf

nind=`bcftools query -l /scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.${i}.GT_unrelatedCMHPO.vcf | wc -l`

python /scratch/pae3g/zaps/software/SelectionHapStats/scripts/H12_H2H1.py /scratch/pae3g/zaps/popgen/${i}.H12.csv $nind -o /scratch/pae3g/zaps/popgen/${i}.CMHPO.200.50.H12results.txt -w 200 -j 50 -d 0
