#module load gcc/7.1.0  openmpi/3.1.4 R/3.6.3; R


library(SeqArray)

seqVCF2GDS("/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.daphnid.whatshap.vcf",
          "/scratch/aob2x/daphnia_hwe_sims/popPhase/shapeitOut/MapJune2020_ann.filter.daphnid.whatshap.gds", storage.option="ZIP_RA")


#seqVCF2GDS("/scratch/aob2x/dest/dest.June14_2020.ann.vcf", "/scratch/aob2x/dest.June14_2020.ann.gds", storage.option="ZIP_RA")

args = commandArgs(trailingOnly=TRUE)
vcf.fn=args[[1]]
gds.fn=gsub(".vcf", ".gds", vcf.fn)

vcf.fn=paste(vcf.fn, ".gz", sep="")

seqVCF2GDS(vcf.fn, gds.fn, storage.option="ZIP_RA")
