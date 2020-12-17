# Create sample list for vcfs
# 10.21.2020

# Library
library(data.table)

# name
name<-read.delim('samps.txt')

dt <- as.data.table(expand.grid(name))
setnames(dt, names(dt), c("name"))
dt <- cbind(data.table(slurmID=c(1:dim(dt)[1])), dt)

write.csv(dt, quote=F, row.names=F,
         file="/project/berglandlab/connor/psmc/output/sample_list")
