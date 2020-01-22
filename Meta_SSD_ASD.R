
library(metafor)

# read in data for lower level soc cog
low_meta_data <- read.csv("/projects/loliver/Systematic_Review/Lower_Data_2020-01-09.csv")
#low_meta_data <- low_meta_data[low_meta_data$Outcome.1!="RMET",]
hi_meta_data <- read.csv("/projects/loliver/Systematic_Review/Higher_Data_2020-01-09.csv")
rmet_meta_data <- read.csv("/projects/loliver/Systematic_Review/RMET_Data_2020-01-09.csv")

# calculate individual effect sizes
# specifically, calculate standardized mean effect sizes (Hedge's g) from the means, SDs, and Ns for the SSD and ASD groups for each paper
# vtype="LS" (the default) uses the usual large-sample approximation to compute the sampling variances
# vtype="UB" provides unbiased estimates of the sampling variances

# low without RMET
eff_low <- escalc(measure="SMD", n1i=Task.SSD.N, n2i=Task.ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
       sd2i=O1.ASD.SD, data = low_meta_data, vtype = "LS", append = TRUE, var.names=c("eff_size","var"))

# effect sizes for moderators of interest (age, IQ, education)
eff_low$age_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.Age.Mean, m2i=ASD.Age.Mean, sd1i=SSD.Age.SD, 
                  sd2i=ASD.Age.SD, data = eff_low, vtype = "LS", append = FALSE, var.names=c("age_diff","var"))[1]
eff_low$IQ_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.IQ.Mean, m2i=ASD.IQ.Mean, sd1i=SSD.IQ.SD, 
                    sd2i=ASD.IQ.SD, data = eff_low, vtype = "LS", append = FALSE, var.names=c("IQ_diff","var"))[1]
eff_low$ed_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.Ed.Mean, m2i=ASD.Ed.Mean, sd1i=SSD.Ed.SD, 
                  sd2i=ASD.Ed.SD, data = eff_low, vtype = "LS", append = FALSE, var.names=c("ed_diff","var"))[1]

# diff in prop male and prop on antipsychotics variables
eff_low$Diff.on.Antipsych <- (eff_low$SSD.Prop.on.Antipsych)-(eff_low$ASD.Prop.on.Antipsych)
eff_low$Diff.Male <- (eff_low$SSD.Prop.Male)-(eff_low$ASD.Prop.Male)


# fit model (also sig with RMET studies included, including with perm testing)
# Restricted maximum-likelihood (REML) estimation is used by default when estimating τ2 
# the REML estimator is approximately unbiased and quite efficient; see Viechtbauer 2005; Veroniki et al., 2015)
# DL (adequate and most commonly used) and PM (most recommended by Veroniki et al.) estimators produce extremely similar results 
fit_low <- rma(yi=eff_size, vi=var, data=eff_low, method="REML")

# with pub year as moderator - sig (accounts for 70% of heterogeneity)
fit_low_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_low)

# with QA Score and age as moderators (more exploratory) 
# did check age and QA alone and NS, as with together
fit_low_contmods <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_low)

# exploratory - prop on antipsychotics - NS
fit_low_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_low)

# exploratory - IQ (full-scale or WRAT) k = 4 (don't use)
fit_low_IQ <- rma(yi=eff_size, vi=var, mods=IQ_diff, data=eff_low[eff_low$IQ.Type=="Full-scale"|eff_low$IQ.Type=="WRAT",])

# exploratory - prop male
fit_low_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_low)

# sensitivity analysis - data availability - NS but use
fit_low_avail <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_low)


#setseed=999
#permutation testing (non-normality of observed effects)
#fit_low_p <- permutest(fit_low, exact=T) # still sig

# get confidence intervals
confint(fit_low)

# forest plot (effects)
forest(fit_low, slab = paste(low_meta_data$Author, low_meta_data$Year, sep = ", "),
         ilab = cbind(low_meta_data$Task.SSD.N, low_meta_data$Task.ASD.N), ilab.xpos = c(-9.5, -6), cex = 0.75)

# funnel and radial plots (pub bias)
funnel(fit_low, main = "Random-Effects Model")

# other plot options for detecting bias
radial(fit_low, main = "Random-Effects Model")
qqnorm(fit_low, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry (pub bias)
# One may be able to detect such asymmetry by testing whether the observed outcomes (or residuals from a model with moderators) 
# are related to their corresponding sampling variances, standard errors (the default, "sei"), or more simply, sample sizes
regtest(fit_low, model="rma", predictor = "sei")

# trim and fill method 
# nonparametric (rank-based) data augmentation technique used to estimate the number of studies missing from a meta-analysis
trimfill(fit_low) # if sig, can generate funnel plot of this model object

# check out influence of specific studies
inf_low <- influence(fit_low)
plot(inf_low, plotdfb = TRUE)

# can also run leave one out analysis - stil sig with Bolte & Poustka out
leave1out(fit_low)


# higher-level soc cog
eff_hi <- escalc(measure="SMD", n1i=Task.SSD.N, n2i=Task.ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
                  sd2i=O1.ASD.SD, data = hi_meta_data, vtype = "LS", append = TRUE, var.names=c("eff_size","var"))

# effect sizes for moderators of interest (age, IQ, education)
eff_hi$age_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.Age.Mean, m2i=ASD.Age.Mean, sd1i=SSD.Age.SD, 
                           sd2i=ASD.Age.SD, data = eff_hi, vtype = "LS", append = FALSE, var.names=c("age_diff","var"))[1]
eff_hi$IQ_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.IQ.Mean, m2i=ASD.IQ.Mean, sd1i=SSD.IQ.SD, 
                          sd2i=ASD.IQ.SD, data = eff_hi, vtype = "LS", append = FALSE, var.names=c("IQ_diff","var"))[1]
eff_hi$ed_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.Ed.Mean, m2i=ASD.Ed.Mean, sd1i=SSD.Ed.SD, 
                          sd2i=ASD.Ed.SD, data = eff_hi, vtype = "LS", append = FALSE, var.names=c("ed_diff","var"))[1]

# diff in prop male and prop on antipsychotics variables
eff_hi$Diff.on.Antipsych <- (eff_hi$SSD.Prop.on.Antipsych)-(eff_hi$ASD.Prop.on.Antipsych)
eff_hi$Diff.Male <- (eff_hi$SSD.Prop.Male)-(eff_hi$ASD.Prop.Male)


# model fit
fit_hi <- rma(yi=eff_size, vi=var, data = eff_hi)

# without Waris et al.
eff_hi_15 <- eff_hi[eff_hi$Author!="Waris et al.",]
fit_hi_15 <- rma(yi=eff_size, vi=var, data = eff_hi_15)

# with pub year as moderator # no sig moderators (w eff_hi or eff_hi_15) except
fit_hi_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_hi)

# with QA Score and age as moderators
fit_hi_contmods <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_hi)

# exploratory - prop on antipsychotics - sig, k = 8; Waris not included either way
fit_hi_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_hi)

# exploratory - IQ (full-scale or WRAT); Waris not included either way
fit_hi_IQ <- rma(yi=eff_size, vi=var, mods=IQ_diff, data=eff_hi[eff_hi$IQ.Type=="Full-scale"|eff_hi$IQ.Type=="WRAT",])

# exploratory - prop male - sig without Waris
fit_hi_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_hi_15)

# stim type - NS; Waris not included either way
fit_hi_stim <- rma(yi=eff_size, vi=var, mods=~factor(Stim.Type), data=eff_hi)

# with verb and vis only (no both) - NS; Waris not included either way
fit_hi_stim_2 <- rma(yi=eff_size, vi=var, mods=~factor(Stim.Type), data=eff_hi[eff_hi$Stim.Type!="both",])

# sensitivity analysis - data availability - NS but use
fit_hi_avail <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_hi)


# permutation testing
#setseed=999
#fit_hi_p <- permutest(fit_hi, exact=F) # still non-sig

# get confidence intervals
confint(fit_hi)

# forest plot (effects)
forest(fit_hi, slab = paste(hi_meta_data$Author, hi_meta_data$Year, sep = ", "),
       ilab = cbind(hi_meta_data$Task.SSD.N, hi_meta_data$Task.ASD.N), ilab.xpos = c(-9.5, -6), cex = 0.75)

# funnel and radial plots (pub bias)
funnel(fit_hi, main = "Random-Effects Model")

# other plot options for detecting bias
radial(fit_hi, main = "Random-Effects Model")
qqnorm(fit_hi, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry (pub bias)
regtest(fit_hi, model="rma", predictor = "sei")

# trim and fill method 
trimfill(fit_hi)

# check out influence of specific studies
inf_hi <- influence(fit_hi)
plot(inf_hi, plotdfb = TRUE)

# leave one out analysis
leave1out(fit_hi)


# RMET
eff_rmet <- escalc(measure="SMD", n1i=Task.SSD.N, n2i=Task.ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
                  sd2i=O1.ASD.SD, data = rmet_meta_data, vtype = "LS", append = TRUE, var.names=c("eff_size","var"))

# effect sizes for moderators of interest (age, IQ, education)
eff_rmet$age_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.Age.Mean, m2i=ASD.Age.Mean, sd1i=SSD.Age.SD, 
                          sd2i=ASD.Age.SD, data = eff_rmet, vtype = "LS", append = FALSE, var.names=c("age_diff","var"))[1]
eff_rmet$IQ_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.IQ.Mean, m2i=ASD.IQ.Mean, sd1i=SSD.IQ.SD, 
                         sd2i=ASD.IQ.SD, data = eff_rmet, vtype = "LS", append = FALSE, var.names=c("IQ_diff","var"))[1]
eff_rmet$ed_diff <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=SSD.Ed.Mean, m2i=ASD.Ed.Mean, sd1i=SSD.Ed.SD, 
                         sd2i=ASD.Ed.SD, data = eff_rmet, vtype = "LS", append = FALSE, var.names=c("ed_diff","var"))[1]

# diff in prop male and prop on antipsychotics variables
eff_rmet$Diff.on.Antipsych <- (eff_rmet$SSD.Prop.on.Antipsych)-(eff_rmet$ASD.Prop.on.Antipsych)
eff_rmet$Diff.Male <- (eff_rmet$SSD.Prop.Male)-(eff_rmet$ASD.Prop.Male)

# model fit
fit_rmet <- rma(yi=eff_size, vi=var, data = eff_rmet)

# with pub year as moderator
fit_rmet_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_rmet)

# with QA Score and age as moderators
fit_rmet_contmods <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_rmet)

# exploratory - prop on antipsychotics
fit_rmet_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_rmet)

# exploratory - IQ (full-scale or WRAT)
fit_rmet_IQ <- rma(yi=eff_size, vi=var, mods=IQ_diff, data=eff_rmet[eff_rmet$IQ.Type=="Full-scale"|eff_rmet$IQ.Type=="WRAT",])

# exploratory - prop male
fit_rmet_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_rmet)

# sensitivity analysis - data availability - NS but use
fit_rmet_avail <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_rmet)


# permutation testing
#setseed=999
#fit_rmet_p <- permutest(fit_rmet, exact=F) # still non-sig

# get confidence intervals
confint(fit_rmet)

# forest plot (effects)
forest(fit_rmet, slab = paste(rmet_meta_data$Author, rmet_meta_data$Year, sep = ", "),
       ilab = cbind(rmet_meta_data$Task.SSD.N, rmet_meta_data$Task.ASD.N), ilab.xpos = c(-9.5, -6), cex = 0.75)

# funnel and radial plots (pub bias)
funnel(fit_rmet, main = "Random-Effects Model")

# other plot options for detecting bias
radial(fit_rmet, main = "Random-Effects Model")
qqnorm(fit_rmet, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry (pub bias)
regtest(fit_rmet, model="rma", predictor = "sei")

# trim and fill method 
trimfill(fit_rmet)

# check out influence of specific studies
inf_rmet <- influence(fit_rmet)
plot(inf_rmet, plotdfb = TRUE)

# leave one out analysis
leave1out(fit_rmet)


