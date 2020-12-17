# Compile the SMC++ output 
# 11/14/2020 

# Libraries
library(foreach)
library(data.table)
library(tidyverse)
library(rjson)
library(igraph)
library(viridis)
library(igraph)

# SMC++ Output
setwd('C:/Users/Conno/Desktop/Fall2020/Berglandlab/smcpp')

# Read in psmc output
data <- data.table(read.csv(file = "split4.csv"))
data$comp <- data[,list(comp=paste(clone1, clone2, sep = " ")),]

# Metadata sample file
samps <- fread("C:/Users/Conno/Desktop/Fall2020/Berglandlab/Data/cloneinfo.txt")
samps <- data.table(samps %>% filter(clone%in%data$samp1 & clone%in%data$samp2))

# Read in psmc output
list <- data.table(read.csv(file = "joint_split.csv"))
list$comp <- list[,list(comp=paste(SC.x, SC.y, sep = " ")),]

# Quality control
#write.csv(list[which(list$comp %in% data$comp == FALSE)], file = "remaining_splits", quote = F, row.names = F)
data$comp[which(data$comp %in% list$comp == FALSE)]

# Make symmetric
g <- graph.data.frame(data, directed=FALSE)
g <- get.adjacency(g, attr="split.year", sparse=FALSE)
g <- data.table(melt(g, value.name = "split.year"))
g[split.year==0]$split.year <- "NA"

# Merge with metadata
new <- data.table(g %>% 
                    left_join(samps %>% select(SC, clone, population, year), by=c("Var1"="SC")) %>% 
                    left_join(samps %>% select(SC, clone, population, year), by=c("Var2"="SC")))

new$year <- paste(new$year.x, new$year.y)
new$pond <- paste(new$population.x, new$population.y)

new <- new[!Var1%in%c("R","E") & !Var2 %in% c("R", "E")]

# Split
ggplot(new, aes(x=Var1, y=Var2, fill=split.year)) +
  geom_tile() +
  labs(x="Clonal lineage", y="Clonal lineage", fill="Split year") +
  theme_classic() + 
  scale_fill_viridis(option = "magma") +
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) 

# Pond
ggplot(new, aes(x=Var1, y=Var2, fill=split.year)) +
  geom_tile() +
  facet_wrap(~pond, scales = "free") +
  labs(x="Clonal lineage", y="Clonal lineage", fill="Split year") +
  theme_classic() + 
  scale_fill_viridis(option = "magma") +
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) 

# Split pond
ggplot(na.omit(new), aes(x=population.x, y=population.y, fill=split.year)) +
  geom_raster() +
  labs(x="Pond", y="Pond", fill="Split year") +
  theme_classic() + 
  scale_fill_viridis(option = "magma") +
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) 
