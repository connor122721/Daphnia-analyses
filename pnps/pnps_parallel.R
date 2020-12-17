# 11.10.2020

module load intel/18.0 intelmpi/18.0 R/3.6.3; R

# Libraries
suppressMessages(library(data.table))
suppressMessages(library(SeqArray))
suppressMessages(library(foreach))
suppressMessages(library(tidyverse))

# set working directory
setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/June2020/")

# Executable in command line
arg <- commandArgs(TRUE)

# Execute in parallel
SC.names <- as.character(arg[1])
samp.n <- as.character(arg[2])
out.name <- as.character(arg[3])
rep <- as.character(arg[4])
tmpdir <- as.character(arg[5])

SC.names="OO"; samp.n=8; out.name="test"; rep=23; tmpdir="test"; year=17

# Load meta-data file
samps <- fread("Superclones201617182019withObtusaandPulicaria_kingcorr_20200623_wmedrd.txt")
samps <- samps[Nonindependent==0][is.na(LowReadDepth)]

# Conditional OO control
if(year == 17) {
samps <- samps[SC%in%c("OO")][year%in%c(2017)][population=="DBunk"]
} else
samps <- samps[SC%in%c("OO")][year%in%c(2018, 2019)][population=="D8"]

# Load GDS file
genofile <- seqOpen("MapJune2020_ann.seq.gds")

# Good snps and chromosomes
snpFilter <- fread("finalsetsnpset01pulex_table_20200623")

# Allele frequency pnps function
pnps.fun <- function(SC.i, sample.n, rep=rep) {

    seqResetFilter(genofile)
    seqSetFilter(genofile, sample.id=sample(samps[SC==SC.i]$clone, sample.n),
               variant.id=snpFilter$variant.ids)

    message("Getting allele counts")
    snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                         alleleCount=seqAlleleCount(genofile, ref.allele=1L))
    
    # Annotations
    message("Annotating variants")
    seqSetFilter(genofile, variant.id=snp.dt$variant.id)

    # Gets annotations
    tmp <- seqGetData(genofile, "annotation/info/ANN") 
    
    # Gets allele frequencies
    tmp.af <- data.table(variant.id= snp.dt$variant.id,
                         af= seqAlleleFreq(genofile, ref.allele=1L))
    
    len1 <- tmp$length
    len2 <- tmp$data

    snp.dt1 <- data.table(len=rep(len1, times=len1),
                          ann=len2,
                          id=rep(snp.dt$variant.id, times=len1))

    # Extracting data between the 2nd and third | symbol
    snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]

    # Collapsing additional annotations to original SNP vector length
    snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=",")),
                           list(variant.id=id)]
    
    snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
    
    # Merges allele frequencies
    snp.dt1.an <- merge(snp.dt1.an, tmp.af, by="variant.id")

    # Final merge
    m <- merge(snp.dt, snp.dt1.an, by="variant.id")

    # Summarize
    m.ag <- m[col%in%c("missense_variant", "synonymous_variant"),
                                     list(pn = sum(col=="missense_variant"), 
                                          ps = sum(col=="synonymous_variant"), 
                                          SC = SC.i, rep = rep, sample.n=sample.n, 
                                          year = year), 
                                     list(af = af,  MAC = alleleCount)]
    
    # Return
    return(m.ag)
}

# Run pnps.out function
pnps.out <- pnps.fun(SC.i=SC.names, sample.n=samp.n, rep=rep)

# Finish message
message(paste("SC", SC.names, "Rep #", rep, sep=" "))

# Saving pnps output
write.csv(pnps.out, file=paste(tmpdir, out.name, SC.names, rep, ".csv", sep=""), quote = F, row.names = F)