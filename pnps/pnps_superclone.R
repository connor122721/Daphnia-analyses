# pnps ratios across clonal lineages
# Connor Murray
# 5.17.2020

# Libraries
library(data.table)
library(SeqArray)
library(foreach)
library(tidyverse)

# set working directory
setwd("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/")

# Load meta-data file
samps <- fread("/project/berglandlab/Karen/MappingDec2019/WithPulicaria/Superclones201617182019withObtusaandPulicaria_kingcorr_20200402_wmedrd.txt")

# Load GDS file
genofile <- seqOpen("MapDec19PulexandObtusaandPulicaria_filtsnps10bpindels_snps_filter_pass_lowGQmiss_ann.seq.gds")

# SNP filter file
snpFilter <- fread("finalsetsnpset01pulex_table_wpulicaria_20200401")

# Allele frequency pnps function
pnps.fun <- function(SC.i, MAC, sample.n) {
    # SC.i="A"; MAC=1; sample.n=25; j=1

    seqResetFilter(genofile)
    seqSetFilter(genofile, 
                 sample.id=sample(samps[population%in%c("D8", "DBunk")][SC==SC.i][Nonindependent==0][year==2017]$clone, sample.n), 
                 variant.id=snpFilter$variant.ids)

    message("Getting allele counts")
    snp.dt <- data.table(variant.id=seqGetData(genofile, "variant.id"),
                         alleleCount=seqAlleleCount(genofile, ref.allele=1L))

    # Annotations
    message("Annotating variants")

    snp.dt <- snp.dt[alleleCount==MAC]
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
    m.ag <- m[alleleCount==MAC][,list(pn=sum(col=="missense_variant"), 
                                      ps=sum(col=="synonymous_variant"), 
                                      SC=SC.i, rep=j), 
                                 list(MAC=alleleCount, af=af)]
    
    # Return
    return(m.ag)
}

# Clones of interest
clones <- unique(samps[population %in% c("D8", "DBunk") & year == "2017" & SC %in% c("A", "C", "E", "F", "OO") & Nonindependent==0]$SC)

# Run pnps.out function
pnps.out <- foreach(j = 1:100, .combine="rbind") %do% {
  o <- foreach(i = clones, .combine="rbind") %do% {
                 message(paste("Sample", i, "#", j, sep=" "))
                 pnps.fun(SC.i=i, MAC=1, sample.n=15)
               }
}

# Saving pnps output
save(pnps.out, file="/project/berglandlab/connor/pnps/pnps.15samp.100rep.2017yr.Rdata")

### Graphing & Analyses ###

# pnps Rdata
load("C:/Users/Conno/Desktop/Spring2020/pnps.daphnia/data/pnps.15samp.100rep.2017yr.Rdata")

# adds pnps and filters 
pnps.out[,pnps:=pn/ps]
pnps.out <- na.omit(pnps.out)
pnps.out <- pnps.out[!pnps > 100]
pnps.out.ag <- pnps.out[,list(pnps= mean(pnps), 
                              lci= quantile(pnps, .025), uci= quantile(pnps, .975),
                              var = var(pnps)), 
                         list(SC,af)]

# pnps at several mutation classes boxplot
ggplot(data=pnps.out[af < 0.042 & SC %in% c("A", "C", "OO")], aes(x= as.factor(af), y= pnps, fill=SC)) +
  geom_boxplot() +
    scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N", "4/2N")) +
    scale_fill_discrete(name = "Clonal lineage", labels = c("A", "C", "Outbred")) +
    labs(x="Mutation frequencies", y="pN/pS") +
    theme_bw() 

# pnps at several mutation classes boxplot across all clones
ggplot(data=pnps.out[af < 0.04], aes(x= as.factor(af), y= pnps, fill=SC)) +
  geom_boxplot() +
  scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N")) +
  scale_fill_discrete(name = "Clonal lineage") +
  labs(x="Mutation frequencies", y="pN/pS", title="2017 DBunk and D8 clones") +
  theme_bw() 

# A vs OO pnps
ggplot(data=pnps.out[af < 0.042 & SC %in% c("A", "OO")], 
       aes(x= as.factor(af), y= pnps, color=SC, group=SC)) +
  geom_jitter(alpha=0.3, width = 0.25) +
  geom_smooth(method = "glm", position = "identity", se = FALSE) +
  scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N", "4/2N")) +
  scale_color_discrete(name = "Clonal lineage", labels = c("A", "Outbred")) +
  labs(x="Mutation frequencies", y="pN/pS", title="2017 DBunk and D8 clones") +
  theme_bw() 

# Facet wrap of all pnps points
ggplot(data=pnps.out[af < 0.042], 
       aes(x= as.factor(af), y= pnps, color=SC, group=SC)) +
  geom_jitter(alpha=0.3, width = 0.2) +
  facet_wrap(~SC, scales = "free") +
  geom_smooth(method = "glm", position = "identity", se = FALSE) +
  scale_x_discrete(labels = c("1/2N", "2/2N", "3/2N", "4/2N")) +
  scale_color_discrete(name = "Clonal lineage", labels = c("A", "C", "E", "F", "Outbred")) +
  labs(x="Mutation frequencies", y="pN/pS") +
  theme_bw() 

# pnps at several mutation classes boxplot 1/2N
ggplot(data=pnps.out[af ==min(af) & SC %in% c("A", "C", "OO")], aes(x= SC, y= pnps, fill=SC)) +
  geom_boxplot()  +
  scale_x_discrete(labels = c("A", "C", "Outbred")) +
  labs(x="Clonal lineage", y="pN/pS") +
  theme_bw() +
  ggtitle("1/2N pN/pS") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

# pnps for 1/2N mutations
ggplot(data=pnps.out.ag[af==min(af)]) +
  geom_point(aes(x=SC, y=pnps, color=SC)) +
  geom_errorbar(aes(x=SC, ymin=lci, ymax=uci, color=SC), width=0.3) +
    labs(x="Clonal lineage", y="pN/pS", title= "1/2N Mutation freuquency") +
    theme_classic() +
    theme(legend.position="none")

