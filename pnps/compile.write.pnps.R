# Compile and output pnps data 
# 11/11/2020 

# Libraries
library(foreach)
library(data.table)
library(tidyverse)

# pnps Output
setwd('/project/berglandlab/connor/pnps/data/output')

# Extract all pnps files
tbl <- list.files(pattern = "*.csv$")

# Read in files as table
out <- foreach(i=1:length(tbl), .combine="rbind") %do% {
  
  # Read file
  X <- read.csv(file=tbl[i], header = T, stringsAsFactors = F)
  
  # Finish
  return(X)
  
}

write.csv(x = out, file = "/project/berglandlab/connor/pnps/data/pnps.out.csv", 
          row.names = F, quote = F)
