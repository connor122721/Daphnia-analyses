# Create parameter files to test slimulations
# 9.8.2020

# Library
library(data.table)

# Exponential start and constant growth population
EG = c(1.5) # exponential growth magnitude
K = c(seq(100, 900, by=100), seq(1000, 9000, by=1000), seq(10000, 100000, by=10000)) # carrying capacity
Rep = 10
nSamp = 15
Gen = c(seq(10, 1000, by=5))

dt <- as.data.table(expand.grid(EG, K, Rep, nSamp, Gen))
setnames(dt, names(dt), c("EG", "K", "Rep", "nSamp", "Gen"))
dt <- cbind(data.table(slurmID=c(1:dim(dt)[1])), dt,
	    simID=round(abs(rnorm(1, 1e8, n=1:dim(dt)[1])), digits=0))

write.csv(dt, quote=F, row.names=F,
         file="/scratch/csm6hg/slim/constant_model_paramList")
