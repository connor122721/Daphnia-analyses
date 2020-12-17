# Compile the SMC++ output 
# 11/20/2020 

module load intel/18.0 intelmpi/18.0 R/3.6.3; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)
library(rjson)

# SMC++ Output
setwd('/scratch/csm6hg/smcpp/split/split')

# Extract all model files
tbl <- list.files(pattern = "*.csv$", recursive=TRUE)

# Read in files as table
out <- foreach(i=1:length(tbl), .combine="rbind") %do% {
  
  # Read in files as table
  out <- data.table(read.table(tbl[i], header=T, sep = ","), 
                    clone1=tstrsplit(tbl[[i]], "_")[[1]],
                    clone2=gsub(tstrsplit(tbl[[i]], "_")[[2]], pattern = ".csv", replacement = ""),
                    #chrom=tstrsplit(tbl[[i]], "\\.")[[2]],
                    csv=tbl[[i]],
                    iteration=i)
  print(i)
  
  # Finish  
  return(out)
  
}

write.csv(x = out, file = "/project/berglandlab/connor/smcpp/smcpp.split4.csv", 
          row.names = F, quote = F)

# Extract all model files
tbl <- list.files(pattern = "*.json$", recursive=TRUE)

# Read in files as table
out <- foreach(i=1:length(tbl), .combine="rbind") %do% {
  
  # Read in files as table
  out <- data.table(readLines(tbl[i]))
  
  split = as.numeric(gsub("[^0-9.-]", "", out[grep(out$V1, pattern = "split")]))
  N0 = as.numeric(gsub("[^0-9.-]", "", out[grep(out$V1, pattern = "N0")][1]))
           
  # Split time  
  if (is.na(split)==TRUE) {
    
    split = as.numeric(gsub("-", "E-", gsub("[^0-9.-]", "", out[grep(out$V1, pattern = "split")])))
    
  }
  
  joint <- data.table(read.csv("/scratch/csm6hg/smcpp/joint_split", header = T))
  
  # Output
  out <- data.table(clone1=tstrsplit(tbl[[i]], "_")[[1]], 
                    clone2=gsub(tstrsplit(tbl[[i]], "_")[[2]], pattern = ".csv", replacement = ""),
                    split=split,
                    N0=N0,
                    split.year=N0*2*split*0.2,
                    split.gen=N0*2*split,
                    csv=tbl[[i]],
                    iteration=i)
  
  p1 <- as.character(joint[which(out[iteration==i]$clone1 == joint$SC.x & out[iteration==i]$clone2 == joint$SC.y)]$samp1)
  p2 <- as.character(joint[which(out[iteration==i]$clone1 == joint$SC.x & out[iteration==i]$clone2 == joint$SC.y)]$samp2)
  
  out$samp1 <- p1
  out$samp2 <- p2
  
  print(i)
  
  # Finish  
  return(out)
  
}

write.csv(x = out, file = "/project/berglandlab/connor/smcpp/split4.csv", 
          row.names = F, quote = F)

# Read in psmc output
data <- data.table(na.omit(read.csv(file = "C:/Users/Conno/Desktop/Fall2020/Berglandlab/smcpp/smcpp.out5.csv")))

# Metadata sample file
samps <- fread("C:/Users/Conno/Desktop/Fall2020/Berglandlab/Data/cloneinfo.txt")
samps <- data.table(samps %>% filter(Nonindependent==0 & is.na(LowReadDepth) & Species=="pulex" & medrd > 12))

# Merge metadata info
data <- data.table(merge(data, samps %>% select(SC, clone, population, year), by = "clone"), data="smc++")

# Read in psmc output
ps <- data.table(na.omit(read.csv(file = "C:/Users/Conno/Desktop/Fall2020/Berglandlab/Data/psmc.vcf1.csv")))

# Merge metadata info
ps$clone <- str_replace_all(ps$name, pattern=".psmc", replacement="")
ps <- data.table(merge(ps, samps %>% select(SC, clone, population, year), by = "clone"), data="psmc/Bed")

ps2 <- data.table(na.omit(read.csv(file = "C:/Users/Conno/Desktop/Fall2020/Berglandlab/Data/psmc.vcf2.csv")))

# Merge metadata info
ps2$clone <- str_replace_all(ps2$name, pattern=".psmc", replacement="")
ps2 <- data.table(merge(ps2, samps %>% select(SC, clone, population, year), by = "clone"), data="psmc")

# SMC Split
split <- data.table(na.omit(read.csv(file = "C:/Users/Conno/Desktop/Fall2020/Berglandlab/smcpp/smcpp.split.csv")))

# Merge metadata info
split <- data.table(merge(split, samps %>% select(SC, clone, population, year), 
                          by.x = "label", by.y="clone"), data="smc++")

# Merge SMC++ & PSMC
tot <- data.table(rbind(ps %>% select(clone, years, eff, data, SC, year, population), 
                        ps2 %>% select(clone, years, eff, data, SC, year, population), 
                        data  %>% select(clone, years, eff, data, SC, year, population)))

# Summary output
ggplot(tot[years>1000 & years<1e5], aes(x=years, y=eff, color=data, group=interaction(clone, data))) +
  # stat_smooth(geom='line', alpha=0.5, se=FALSE, span = 0.4) +
  geom_step(alpha=0.5) +
  facet_wrap(~SC) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(scaled = TRUE) +
  labs(x="Years (log10)", y="Ne (log10)") +
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) 

# Split
ggplot(split[x>1000 & x<1e5][iteration==2], aes(x=x, y=y, group=interaction(label,iteration,as.factor(plot_num)))) +
  # stat_smooth(geom='line', alpha=0.5, se=FALSE, span = 0.4) +
  geom_step(alpha=0.5) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks(scaled = TRUE) +
  labs(x="Years (log10)", y="Ne (log10)") +
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18)) 
