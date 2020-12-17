# Libraries
library(data.table)
library(foreach)
library(tidyverse)
library(cowplot)
library(ggforce)
library(ggridges)

# Working directory
setwd("C:/Users/Conno/Desktop/Fall2020/Berglandlab/V-stats")

# Load data
load(file = "phenos_F1.Rdata")

# Epphipia data
ggplot(data=epp[SC %in% c("A", "C", "AxCF1", "selfedC", "selfedA")], 
       aes(x=fill, y=SC, height=..density.., fill=SC, group=SC)) +
   geom_density_ridges(alpha=0.7) +
   labs(x="Proportion ephippia filled", y="") +
   theme_classic() + 
   theme(axis.text.x = element_text(face="bold", size=16),
         axis.text.y = element_text(face="bold", size=16),
         axis.title.x = element_text(face="bold", size=18),
         axis.title.y = element_text(face="bold", size=18),
         legend.background = element_blank(),
         legend.text = element_text(size=14),
         legend.title = element_text(face="bold", size=16),
         legend.position = "none") 

# Male data
ggplot(data=male[SC %in% c("A", "C", "AxCF1", "selfedC", "selfedA")], 
       aes(x=propmale, y=SC, height=..density.., fill=SC, group=SC)) +
   geom_density_ridges(alpha=0.7) +
   labs(x="Proportion males", y="") +
   theme_classic() + 
   theme(axis.text.x = element_text(face="bold", size=16),
         axis.text.y = element_text(face="bold", size=16),
         axis.title.x = element_text(face="bold", size=18),
         axis.title.y = element_text(face="bold", size=18),
         legend.background = element_blank(),
         legend.text = element_text(size=14),
         legend.title = element_text(face="bold", size=16),
         legend.position = "none") 

# V statistic for ephippia

numer = var(par) - (var(par1)/4*n(par1)) - (var(par2)/4*n(par2))

denom = var(F2)*H2*c

v = numer/denom








nuc.div <- na.omit(sim[,list(avg.pi = mean(tP.stan),
                                avg.watterson = mean(tW.stan),
                                avg.tajima = mean(tajimaD)),
                           list(EG, K, gen)])

sort(unique(sim$gen))
sort(unique(sim$K))

# Labels for heatmaps
label <- data.table(x=c(38, 9, 47, 79),
                    y=c(14, 15, 14, 14),
                    label=c("A", "B", "C", "D"))

# Theta pi
a <- ggplot() +
   geom_tile(data=nuc.div[gen < 505], 
             aes(x=as.factor(gen), y=as.factor(K), fill=avg.pi)) +
   geom_ellipse(aes(x0=38, y0=14, a=7, b=4, angle=0, fill=2.1e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=9, y0=15, a=7, b=6, angle=0, fill=5.3e-7), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=47, y0=14, a=9, b=3, angle=0, fill=2.5e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=79, y0=14, a=10, b=3, angle=0, fill=4.4e-6), 
                size=1, color="white", alpha=0.8) +
   geom_label(data=label, aes(x=x,y=y,label=label)) +
   labs(x="Exponential growth", y="K", fill="Theta pi") +
   scale_fill_viridis(option = "A") +
   scale_y_discrete(breaks=c(100, 500, 1000, 5000, 10000, 50000, 80000)) +
   theme_classic() + 
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_blank(),
         axis.title.y = element_text(face="bold", size=10),
         axis.ticks.x = element_blank()) 

# Theta Watterson
b <- ggplot() +
   geom_tile(data=nuc.div[gen < 505], 
             aes(x=as.factor(gen), y=as.factor(K), fill=avg.watterson)) +
   geom_ellipse(aes(x0=38, y0=14, a=7, b=4, angle=0, fill=7.1e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=9, y0=15, a=7, b=6, angle=0, fill=1.7e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=47, y0=14, a=9, b=3, angle=0, fill=8.3e-6), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=79, y0=14, a=10, b=3, angle=0, fill=1.3e-5), 
                size=1, color="white", alpha=0.8) +
   geom_label(data=label, aes(x=x,y=y,label=label)) +
   labs(x="Exponential growth", y="K", fill="Theta Wa") +
   scale_fill_viridis(option = "A") +
   scale_y_discrete(breaks=c(100, 500, 1000, 5000, 10000, 50000, 80000)) +
   theme_classic() + 
   theme(axis.text.x = element_blank(),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_blank(),
         axis.title.y = element_text(face="bold", size=10),
         axis.ticks.x = element_blank())

# Tajima's D
c <- ggplot() +
   geom_tile(data=nuc.div[gen < 505], 
             aes(x=as.factor(gen), y=as.factor(K), fill=avg.tajima)) +
   geom_ellipse(aes(x0=38, y0=14, a=7, b=4, angle=0, fill=-2.70), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=9, y0=15, a=7, b=6, angle=0, fill=-2.55), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=47, y0=14, a=9, b=3, angle=0, fill=-2.66), 
                size=1, color="white", alpha=0.8) +
   geom_ellipse(aes(x0=79, y0=14, a=10, b=3, angle=0, fill=-2.59), 
                size=1, color="white", alpha=0.8) +
   geom_label(data=label, aes(x=x,y=y,label=label)) +
   labs(x="Generation", y="K", fill="Tajima's D") +
   scale_fill_viridis(option = "A") +
   scale_x_discrete(breaks=c(10, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500)) +
   scale_y_discrete(breaks=c(100, 500, 1000, 5000, 10000, 50000, 80000)) +
   theme_classic() + 
   theme(axis.text.x = element_text(face="bold", size=8, vjust = 0.4),
         axis.text.y = element_text(face="bold", size=8),
         axis.title.x = element_text(face="bold", size=10),
         axis.title.y = element_text(face="bold", size=10)) 

plot_grid(a,b,c, nrow=3)

