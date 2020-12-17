library(data.table)
library(foreach)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(rehh)

y<-foreach(i=c(1:40), .errorhandling="remove")%do% {
    load(paste0("/scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.scaffold", i, ".replaced_unrelatedCMHPO_ihs.Rdat"))
    return(as.data.table(scan.ihs$ihs))
}

y<-rbindlist(y)
y[,markernum:=c(1:nrow(y))]

ggplot(y, aes(x=markernum, y=IHS, color=CHR))+geom_point()+theme(legend.position = "none")

y<-y[!is.na(IHS)]
y[order(IHS)]

#what does dgrp ihs look like?
dgrp<-fread("/scratch/pae3g/revisions/evolution/ihs.txt")

dgrp.ihs<-ihh2ihs(dgrp)

dgrp.ihs.dt<-as.data.table(dgrp.ihs$ihs)
dgrp.ihs.dt[,markernum:=c(1:nrow(dgrp.ihs.dt))]

ggplot(dgrp.ihs.dt, aes(x=markernum, y=IHS, color=CHR))+geom_point()+theme(legend.position = "none")

#so these numbers are higher than what we see in dgrp

#take a look at scaffold11   103672 

i=11
hh <- data2haplohh(hap_file = paste0("/scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.scaffold", i, ".replaced_unrelatedCMHPO.vcf"),
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")
scan <- scan_hh(hh)

scan.dt<-as.data.table(scan)
scan.dt[,markernum:=c(1:nrow(scan.dt))]
scan.dt[POSITION==103672]
ehh <- calc_ehh(hh, mrk=552 )
plot(ehh)
hh_subset <- subset(hh, select.mrk = 452:652)
pdf("/scratch/pae3g/zaps/ihs_sc11_103672.pdf", height=8, width=8)
plot(
    hh_subset,
    mrk = 101,
    group_by_allele = TRUE,
    ignore.distance = TRUE,
    col = c(NA, "red"),
    linecol = c("lightblue", "lightpink"),
    mrk.col = "black",
    cex = 0.1,
    pos.lab.hap = "none",
    pos.lab.mrk = "none"
)
dev.off()

#scaffold3  1700193
i=3
hh <- data2haplohh(hap_file = paste0("/scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.scaffold", i, ".replaced_unrelatedCMHPO.vcf"),
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")
scan <- scan_hh(hh)
scan.dt<-as.data.table(scan)
scan.dt[,markernum:=c(1:nrow(scan.dt))]
scan.dt[POSITION==1700193]
ehh <- calc_ehh(hh, mrk=55095 )
plot(ehh)

hh_subset <- subset(hh, select.mrk = 54995:55195)
pdf("/scratch/pae3g/zaps/ihs_sc13_1700193.pdf", height=8, width=8)
plot(
    hh_subset,
    mrk = 101,
    group_by_allele = TRUE,
    ignore.distance = TRUE,
    col = c(NA, "red"),
    linecol = c("lightblue", "lightpink"),
    mrk.col = "black",
    cex = 0.1,
    pos.lab.hap = "none",
    pos.lab.mrk = "none"
)
dev.off()

#  scaffold4   481724 

i=4
hh <- data2haplohh(hap_file = paste0("/scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.scaffold", i, ".replaced_unrelatedCMHPO.vcf"),
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

hh.tab<-data.table(chr=chr.name(hh),
                   pos=positions(hh),
                   markernum=c(1:length(positions(hh))))

hh.tab[pos==481724]
ehh <- calc_ehh(hh, mrk=7252 )
plot(ehh)

hh_subset <- subset(hh, select.mrk = 7152:7352)
pdf("/scratch/pae3g/zaps/ihs_sc4_481724.pdf", height=8, width=8)
plot(
    hh_subset,
    mrk = 101,
    group_by_allele = TRUE,
    ignore.distance = TRUE,
    col = c(NA, "red"),
    linecol = c("lightblue", "lightpink"),
    mrk.col = "black",
    cex = 0.1,
    pos.lab.hap = "none",
    pos.lab.mrk = "none"
)
dev.off()


#look at area around Ace
i=1
load(paste0("/scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.scaffold", i, ".replaced_unrelatedCMHPO_ihs.Rdat"))
scan.dt<-as.data.table(scan.ihs$ihs)

scan.dt[POSITION>1905000&POSITION<1942000]
ggplot(scan.dt[POSITION>1905000&POSITION<1942000], aes(x=POSITION, y=IHS))+geom_point()

#ihs vs missing data/snp density
y<-foreach(i=c(1:40), .errorhandling="remove")%do% {
    load(paste0("/scratch/pae3g/zaps/popgen/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.scaffold", i, ".replaced_unrelatedCMHPO_ihs.Rdat"))
    return(as.data.table(scan.ihs$ihs))
}

y<-rbindlist(y)
y[,markernum:=c(1:nrow(y))]

miss<-fread("/scratch/pae3g/zaps/haplotype_calling/zaprionus_2020.biallelic.miss80.minQ50.minDP5.maxDP60.unrelated_CMHPO.lmiss")
setnames(miss, "POS", "POSITION")

y<-merge(y, miss, by=c("CHR", "POSITION"))

z<-foreach(snp=c(1:nrow(y)))%dopar%{
    if(snp%%1000==0){
        print(snp)
    }
    m<-y[snp, markernum]
    chrom<-y[snp,CHR]
    subset=y[CHR==chrom & markernum > (m-100) & markernum < (m+100)]
    return(data.table(markernum=m,
                      avg.missing=mean(subset$F_MISS, na.rm=T),
                      marker_density=201/(max(subset$POSITION)-min(subset$POSITION))))
}
z<-rbindlist(z)

save(z, file="/scratch/pae3g/zaps/popgen/ihs_missing_window.Rdat")
#H12/G12
z<-foreach(i=c(1:40), .errorhandling="remove")%do% {
    
    a<-(fread(paste0("/scratch/pae3g/zaps/popgen/scaffold", i, ".CMHPO.200.50.H12results.txt")))
    a[,scaffold:=i]
    return(a)
    }
z<-rbindlist(z)


setnames(z, c("pos", "start", "end", "nhaps", "hfs", "line", "H1", "H2", "H12", "H2_H1", "H123", "chr"))
z[,markernum:=c(1:nrow(z))]

ggplot(z, aes(x=markernum, y=H12, color=as.factor(chr)))+geom_point()+theme(legend.position = "none")
ggplot(z, aes(x=markernum, y=H2_H1, color=as.factor(chr)))+geom_point()+theme(legend.position = "none")



