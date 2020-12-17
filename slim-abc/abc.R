# Libraries
library(abc)
library(abc.data)
library(data.table)
library(foreach)
library(tidyverse)

# Working directory
setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/slim-abc/Data/")

# Empirical diversity data
data <- data.table(na.omit(read.csv("MapJune2020.recode.finalsnpset.F1Wilds.diversity.csv", 
                                    header = T)))

data <- data.table(data %>% select(tW.stan.pop.1, tP.stan.pop.1, tajima.D.pop.1, SC) %>% 
                     mutate(tP.stan=as.numeric(tP.stan.pop.1),
                            tW.stan=as.numeric(as.character(tW.stan.pop.1)),
                            tajimaD=as.numeric(tajima.D.pop.1),
                            SC=SC))

ggplot(data[!SC=="OO"], aes(x=SC, y=tP.stan, color=as.factor(SC))) + geom_boxplot()

# Mean of clonal lineages
data <- data.table(data %>% group_by(SC) %>% select(tP.stan, tW.stan, tajimaD) %>% 
                     summarise_each(funs = mean))

# Old diversity
data2 <- data.table(readRDS("filtered-10rep-15samp.rds"))
colnames(data2)[1:6] <- c("tW", "tW.stan", "tP", "tP.stan", "tP.eff", "tajimaD")

ggplot(data2[!SC=="OO"], aes(x=SC, y=tP.stan, color=as.factor(SC))) + geom_point()

data2 <- data.table(data2 %>% group_by(SC) %>% select(tP.stan, tW.stan, tajimaD)%>% 
                     summarise_each(funs = mean))

# Simulation diversity data
sim <- data.table(rbind(read.csv("constant.output.10.20.csv")))

colnames(sim)[1:6] <- c("tW", "tW.stan", "tP", "tP.stan", "tP.eff", "tajimaD")
sim[, seed:=tstrsplit(sim$vcf, "_")[10]]

#sim <- data.table(sim %>% select(tP.stan, tW.stan, tajimaD, gen, K) %>% group_by(gen, K) %>% summarise_each(funs=mean))
ggplot(sim[EG < 1.4], aes(x=gen, y=tP.stan, color=as.factor(K))) + geom_boxplot()

# Model information
sim[, model:=paste("constant", EG, K, gen, sep="_")]
models <- as.character(sim$model)

# Sim stats
sim.stat <- data.table(sim %>% select(tP.stan, tW.stan, tajimaD))

# ABC using rejection
cv.rej <- abc(target = data[1,-c(1)],
              param=sim[,c(19, 20, 24)], 
              sumstat=sim.stat,
              tol=0.1, 
              method="rejection")

# Parameter inference
summary(cv.rej)

