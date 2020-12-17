# Create parameter files
# 11.12.2020

# Library
library(data.table)

SC = c("A", "B", "C", "D", "E", "F", "H", "I", "J", "K", "L", "OO")
rep = c(1:100)
sample.n = 8
out.name = "MapJune2020.recode.pnps2"

dt <- as.data.table(expand.grid(SC, rep))
setnames(dt, names(dt), c("SC", "rep"))
dt <- cbind(data.table(id=c(1:dim(dt)[1]), sample.n=sample.n, out.name=out.name), dt)

write.csv(dt, quote=F, row.names=F, file="/project/berglandlab/connor/pnps/scripts/SC_model_paramList")
