# Compile, scale, and output PSMC data 
# 11/13/2020 

module load intel/18.0 intelmpi/18.0 R/3.6.3; R

# Libraries
library(foreach)
library(data.table)
library(tidyverse)

# PSMC Output
setwd('/project/berglandlab/connor/psmc/output/vcf2/')

# Extract all PSMC files
tbl <- list.files(pattern = "*.psmc$")

# Read in files as table
out <- foreach(i=1:length(tbl), .combine="rbind", .errorhandling = "pass") %do% {
  
  i.iteration=25; mu=5.69e-09; s=100; g=1/20
  
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

write.csv(x = out, file = "/project/berglandlab/connor/psmc/output/psmc.vcf2.csv", 
          row.names = F, quote = F)

# Read in psmc output
data <- data.table(na.omit(read.csv(file = "C:/Users/Conno/Desktop/Fall2020/Berglandlab/Data/psmc.vcf1.csv")))

# Metadata sample file
samps <- fread("C:/Users/Conno/Desktop/Fall2020/Berglandlab/Data/cloneinfo.txt")
samps <- data.table(samps %>% filter(Nonindependent==0 & is.na(LowReadDepth) & Species=="pulex"))

# Merge metadata info
data$clone <- str_replace_all(data$name, pattern=".psmc", replacement="")
data <- data.table(merge(data, samps %>% select(SC, clone, population, year), by = "clone"))
unique(data$clone)

# Summary output
ggplot(data[years > 1e3], aes(x=years, y=eff, color=SC, group=clone)) +
  geom_step(alpha=0.4) +
  facet_wrap(~SC) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x="Years (log10)", y="Ne (log10)") +
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=16),
        axis.text.y = element_text(face="bold", size=16),
        axis.title.x = element_text(face="bold", size=18),
        axis.title.y = element_text(face="bold", size=18),
        legend.position = "none") 
