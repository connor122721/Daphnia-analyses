# Calculating selection statistics across phased/imputed VCF
# 12.15.2020

# Libraries
library(data.table)
library(rehh)

setwd("/project/berglandlab/connor/rehh/")

# Command elements
arg <- commandArgs(trailingOnly = T)
i = arg[1]

# Message
print(paste("Chromosome", i, sep=": "))

# Subsetted VCF
hh <- data2haplohh(hap_file = "MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.vcf",
                   chr.name = i,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

# Haplotype scan
scan <- scan_hh(hh, threads = 10)

# Compile data
scan.ihs <- data.table(ihh2ihs(scan))

# Save data
save(scan.ihs, file=paste0("MapJune2020_ann.hyrbid_strategy.3species.whatshap.shapeit.", 
                           fread("goodChrom.txt", header = F)$V1[i], 
                           ".ihs.Rdata"))
