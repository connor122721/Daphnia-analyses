# Coalescent based models testing population structure on measures of nucleotide diversity
# Connor Murray 5.9.2020

# Libraries
library(scrm)
library(data.table)
library(coala)
library(tidyverse)
library(cowplot)
library(foreach)

### Model constants and assumptions ###

# 12 diploid chromosomes, 1Mb long, 100 loci, no recombination, 
# population scaled mutation rate = 4*1000000*5.69e-09; 1000000 = Ne(0)

{
  # Starting parameters
  rep = 50 # number of bootstraps per model
  exp.inc.rt = seq(from = 0.05, to = 5, by = 0.1) # Exp. increase rate
  exp.dec.rt = seq(from = -5, to = -0.05, by = 0.1) # Exp. decrease rate
  inst.inc.rt = seq(from = 0.05, to = 5, by = 0.1) # Instant increase pop. size
  inst.dec.rt = seq(from = 5, to = 0.05, by = -0.1) # Instant decrease pop. size
  sample.n = 10 # number of chromosomes
  loci.n = 50 # number of loci
}

# Simulate all models
model.out <- foreach(i = 1:max(rep), .combine="rbind") %do% {
    
    # i = 1
    
    # Exp increasing at t(0)
    exp.inc.mod <- coal_model(sample_size = sample.n, loci_number = loci.n,
                       loci_length = 1000000, ploidy = 2) +
        feat_growth(rate = exp.inc.rt[i], time = 0) +
        feat_growth(rate = 0, time = 1) +
        feat_mutation(rate = 4*1000000*5.69e-09) +
        feat_recombination(rate = 0) +
        sumstat_nucleotide_div()

    # Mean of every bootstrap
    exp.inc <- data.table(pie=sapply(simulate(exp.inc.mod, nsim = rep), function(x) mean(x$pi)),
                          rep = 1:rep, rate = exp.inc.rt[i], model = "exp.inc")

    # Exp decreasing at t(0)
    exp.dec.mod <- coal_model(sample_size = sample.n, loci_number = loci.n,
                     loci_length = 1000000, ploidy = 2) +
        feat_growth(rate = exp.dec.rt[i], time = 0) +
        feat_growth(rate = 0, time = 1) +
        feat_mutation(rate = 4*1000000*5.69e-09) +
        feat_recombination(rate = 0) +
        sumstat_nucleotide_div()

    # Mean of every bootstrap
    exp.dec <- data.table(pie=sapply(simulate(exp.dec.mod, nsim = rep), function(x) mean(x$pi)),
                          rep = 1:rep, rate = exp.dec.rt[i], model = "exp.dec")

    # Instant increase in pop from t(1)
    inst.inc.mod <- coal_model(sample_size = sample.n, loci_number = loci.n,
                     loci_length = 1000000, ploidy = 2) +
        feat_size_change(1, time = 0) +
        feat_size_change(inst.inc.rt[i], time = 1) +
        feat_mutation(rate = 4*1000000*5.69e-09) +
        feat_recombination(rate = 0) +
        sumstat_nucleotide_div()

    # Mean of every bootstrap
    inst.inc <- data.table(pie=sapply(simulate(inst.inc.mod, nsim = rep), function(x) mean(x$pi)),
                           rep = 1:rep, rate = inst.inc.rt[i], model = "inst.inc")

    # Instant decrease in pop from t(1)
    inst.dec.mod <- coal_model(sample_size = sample.n, loci_number = loci.n,
                      loci_length = 1000000, ploidy = 2) +
        feat_size_change(inst.dec.rt[i], time = 0) +
        feat_size_change(1, time = 1) +
        feat_mutation(rate = 4*1000000*5.69e-09) +
        feat_recombination(rate = 0) +
        sumstat_nucleotide_div()

    # Mean of every bootstrap
    inst.dec <- data.table(pie=sapply(simulate(inst.dec.mod, nsim = rep), function(x) mean(x$pi)),
                           rep = 1:rep, rate = inst.dec.rt[i], model = "inst.dec")

    # Combine final
    total <- rbind(exp.inc, exp.dec, inst.inc, inst.dec)

    # Cool % complete message
    message(paste("Replicate # ", i, " : ", i/max(rep)*100, "% complete", sep=""))
    
return(total)
    
}

# Summarize
tot <- data.table(model.out %>% group_by(model, rate) %>%
                    summarise(mean.pie = mean(pie),
                              lq.pie = quantile(pie, probs = 0.025),
                              uq.pie = quantile(pie, probs = 0.925),
                              med.pie = median(pie)))

# Plotting
p1 <- ggplot(tot[model=="exp.inc"], aes(x=rate, y=med.pie, color=model)) +
  geom_point(alpha=0.8, color="purple")  +
  geom_errorbar(aes(x=rate, ymin=lq.pie, ymax=uq.pie), width=0, alpha=0.8, color="purple") +
  geom_smooth(method = "loess", se = FALSE, size = 1, linetype=2, color="purple") +
  labs(x="Rate of exponential increase", y="Nucleotide diversity (pie)",
       title="Exponential population increase") +
  theme_classic() +
  theme(legend.position="none")

p2 <- ggplot(tot[model=="exp.dec"], aes(x=rate, y=med.pie)) +
  geom_point(alpha=0.8, color="green")  +
  geom_errorbar(aes(x=rate, ymin=lq.pie, ymax=uq.pie), width=0, alpha=0.8, color="green") +
  geom_smooth(method = "loess", se = FALSE, size = 1, linetype=2, color="green") +
  labs(x="Rate of exponential decrease", y="Nucleotide diversity (pie)",
       title="Exponential population decrease") +
  theme_classic() +
  theme(legend.position="none")

p3 <- ggplot(tot[model=="inst.inc"], aes(x=rate, y=med.pie)) +
  geom_point(alpha=0.8, color="red")  +
  geom_errorbar(aes(x=rate, ymin=lq.pie, ymax=uq.pie), width=0, alpha=0.8, color="red") +
  geom_smooth(method = "loess", se = FALSE, size = 1, linetype=2, color="red") +
  labs(x="Rate of population increase", y="Nucleotide diversity (pie)", 
       title="Instant population Increase") +
  theme_classic() +
  theme(legend.position="none")

p4 <- ggplot(tot[model=="inst.dec"], aes(x=rate, y=med.pie)) +
  geom_point(alpha=0.8, color="blue")  +
  geom_errorbar(aes(x=rate, ymin=lq.pie, ymax=uq.pie), width=0, alpha=0.8, color="blue") +
  geom_smooth(method = "loess", se = FALSE, size = 1, linetype=2, color="blue") +
  labs(x="Rate of population decrease", y="Nucleotide diversity (pie)", 
       title="Instant population decrease") +
  theme_classic() +
  theme(legend.position="none")

plot_grid(p1,p2,p3,p4)

