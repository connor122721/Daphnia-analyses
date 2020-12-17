# Create parameter files
# 9.8.2020

# Library
library(data.table)

SC = c("A", "B", "C", "D", "E", "F", "OO")
chrom = c(1:12)
rep = c(1:100)
sample.n = 15
upper.maf = 0.25
lower.maf = 0
out.name = "100rep.15samp.diversity"

dt <- as.data.table(expand.grid(SC, chrom, rep))
setnames(dt, names(dt), c("SC", "chrom", "rep"))
dt <- cbind(data.table(id=c(1:dim(dt)[1]), sample.n=sample.n, upper.maf=upper.maf, lower.maf=lower.maf, out.name=out.name), dt)

write.csv(dt, quote=F, row.names=F, file="/project/berglandlab/connor/PopGenome/scripts/SC_model_paramList")

