# Create parameter files
# 10.22.2020

# Library
library(data.table)

#SC = c("A", "B", "C", "D", "E", "F", "OO")
SC="BCF"
chrom = c(1:12)
#rep = c(1:100)
upper.maf = 0.5
lower.maf = 0
out.name = "MapJune2020.recode.finalsnpset.phasedBCF"

dt <- as.data.table(expand.grid(SC, chrom))
setnames(dt, names(dt), c("SC", "chrom"))
dt <- cbind(data.table(id=c(1:dim(dt)[1]), upper.maf=upper.maf, lower.maf=lower.maf, out.name=out.name), dt)

write.csv(dt, quote=F, row.names=F, file="/project/berglandlab/connor/PopGenome/scripts/SC_model_paramList")

