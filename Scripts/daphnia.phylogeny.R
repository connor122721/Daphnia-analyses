# Libraries
library(ape)
library(data.table)
library(tidyverse)
library(phytools)

pul.tree <- force.ultrametric(rtree(n = 4))
pul.tree$tip.label <- c( "D.obtusa", "D.pulicaria", "D.pulex", "Simocephalus")

root(pul.tree, outgroup = "Simocephalus")
plot(pul.tree, font = 4, cex = 2)     

