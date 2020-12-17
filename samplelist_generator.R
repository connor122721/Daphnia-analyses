# Create sample list for vcfs
# 11.18.2020

module load intel/18.0 intelmpi/18.0 R/3.6.3; R

# Library
library(data.table)
library(tidyverse)
library(gtools)

# Get random individual per SC
name <- data.table(read.delim('clone.info.txt', header = T, row.names = NULL))
name <- name %>% filter(Nonindependent==0 & is.na(LowReadDepth) & Species=="pulex" & medrd >= 15)
rand.sc <- c(as.character(data.table(name %>% filter(population %in% c("D8","DBunk","DCat")) %>% 
                        group_by(SC) %>% sample_n(1))$clone))

# Chromosomes of interest
chrom <- c('Scaffold_1863_HRSCAF_2081', 'Scaffold_1931_HRSCAF_2197', 'Scaffold_2158_HRSCAF_2565', 
           'Scaffold_2217_HRSCAF_2652', 'Scaffold_2373_HRSCAF_2879', 'Scaffold_6786_HRSCAF_7541', 
           'Scaffold_7757_HRSCAF_8726', 'Scaffold_9197_HRSCAF_10753', 'Scaffold_9198_HRSCAF_10754', 
           'Scaffold_9199_HRSCAF_10755', 'Scaffold_9200_HRSCAF_10757', 'Scaffold_9201_HRSCAF_10758')

# Removes duplicates
dt <- data.table(combinations(n = 33, r = 2, v = rand.sc, repeats.allowed = FALSE))

# Expand and create sample list
dt <- expand_grid(dt, chrom)
setnames(dt, old=names(dt), new=(c("samp1", "samp2", "chr")))

# Merge with metadata
dt <- data.table(dt %>% 
                   left_join(name %>% select(clone, SC), by=c("samp1"="clone")) %>% 
                   left_join(name %>% select(clone, SC), by=c("samp2"="clone")))
                 
# Remove within clones
dt <- dt[!samp1==samp2][!SC.x==SC.y]
dt <- data.table(cbind(data.table(slurmID=c(1:dim(dt)[1])), dt))

# Write output
write.csv(dt, quote=F, row.names=F,
         file="/project/berglandlab/connor/smcpp/samples_joint")

# Output for split command
new <- data.table(combinations(n = 33, r = 2, v = rand.sc, repeats.allowed = FALSE))
setnames(new, old=names(new), new=(c("samp1", "samp2")))

# Merge with metadata
new <- data.table(new %>% 
                   left_join(name %>% select(clone, SC), by=c("samp1"="clone")) %>% 
                   left_join(name %>% select(clone, SC), by=c("samp2"="clone")))

# Remove within clones
new <- new[!samp1==samp2]
new <- data.table(cbind(data.table(slurmID=c(1:dim(new)[1])), new))

# Write output
write.csv(new, quote=F, row.names=F,
          file="/project/berglandlab/connor/smcpp/joint_split")

