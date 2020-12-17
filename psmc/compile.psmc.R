# Compile, scale, and output PSMC data 
# 11/13/2020 

# Libraries
library(foreach)
library(data.table)
library(tidyverse)

# Commands 
arg <- commandArgs(TRUE)
path.name <- arg[1] # pathway to output
out.name <- arg[2] # the output file name
mu <- arg[3]
s <- arg[4]
g <- arg[5]

# Working directory
setwd(path.name)

# Extract all PSMC files
tbl <- list.files(pattern = "*.psmc$")

# Read in files as table
out <- foreach(i=1:length(tbl), .combine="rbind") %do% {
  
  # Read file
  X <- scan(file=tbl[i], what="", sep="\n", quiet=TRUE)
  START <- grep("^RD",X)
  END <- grep("^//",X)          
  
  # Last iterate
  X <- X[START[i.iteration]:END[i.iteration]]
  
  # Extract info
  TR <- grep("^TR",X,value=TRUE)
  RS <- grep("^RS",X,value=TRUE)
  
  # Make tmp files and extract info
  tmp.name <- paste(tbl[i], "temp.psmc.result", sep = "")
  
  write(TR, tmp.name)
  theta0 <- as.numeric(read.table(tmp.name)[1,2])
  N0 <- theta0/4/mu/s
  
  write(RS, tmp.name)
  a <- read.table(tmp.name)
  Generation <- as.numeric(2*N0*a[,3])
  Ne <- as.numeric(N0*a[,4])
  
  file.remove(tmp.name)
  
  # Scale
  n.points <- length(Ne)
  YearsAgo <- c(as.numeric(rbind(Generation[-n.points],Generation[-1])),
              Generation[n.points])*g
  Ne <- c(as.numeric(rbind(Ne[-n.points],Ne[-n.points])),
        Ne[n.points])
  
   # Message
  print(paste("Sample #", i, "-", round(i/length(tbl), 3)*100, "%", sep=" "))
  
  # Final Output
  fin <- data.table(years=YearsAgo, eff=Ne, name=tbl[i], mu=mu, s=s, g=g, iteration=i)
 
  # Finish
  return(fin)
  
}

# Write output
write.csv(x = out, file = out.name, row.names = F, quote = F)