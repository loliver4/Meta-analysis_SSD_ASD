
library(metafor)

# read in data for lower level soc cog
meta_data_low <- read.csv("/projects/loliver/Systematic_Review/Approved_Paper_Data_Lower.csv")
meta_data_low <- meta_data_low[meta_data_low$Outcome.1!="RMET",]
meta_data_hi <- read.csv("/projects/loliver/Systematic_Review/Approved_Paper_Data_Higher.csv")
meta_data_rmet <- read.csv("/projects/loliver/Systematic_Review/Approved_Paper_Data_RMET.csv")

# calculate individual effect sizes
# specifically, calculate standardized mean effect sizes (Hedge's g) from the means, SDs, and Ns for the SSD and ASD groups for each paper
# vtype="LS" (the default) uses the usual large-sample approximation to compute the sampling variances
# vtype="UB" provides unbiased estimates of the sampling variances

# low without RMET (also sig with RMET studies included, including with perm testing)
eff_low <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
       sd2i=O1.ASD.SD, data = meta_data_low, vtype = "LS", append = FALSE, var.names=c("eff_size","var"))

# fit model
# Restricted maximum-likelihood (REML) estimation is used by default when estimating τ2 
# the REML estimator is approximately unbiased and quite efficient; see Viechtbauer 2005)
fit_low <- rma(yi=eff_size, vi=var, data=eff_low)

# with pub year as moderator
fit_low_yr <- rma(yi=eff_size, vi=var, mods=meta_data_low$Year, data=eff_low)

#setseed=999
#permutation testing (non-normality of observed effects)
#fit_low_p <- permutest(fit_low, exact=T) # still sig

# get confidence intervals
confint(fit_low)

# forest plot (effects)
forest(fit_low, slab = paste(meta_data_low$Author, meta_data_low$Year, sep = ", "),
         ilab = cbind(meta_data_low$SSD.N, meta_data_low$ASD.N), ilab.xpos = c(-9.5, -6), cex = 0.75)

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

# can also run leave one out analysis
leave1out(fit_low)


# higher-level soc cog
eff_hi <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
                  sd2i=O1.ASD.SD, data = meta_data_hi, vtype = "LS", append = FALSE, var.names=c("eff_size","var"))

fit_hi <- rma(yi=eff_size, vi=var, data = eff_hi)

# permutation testing
#setseed=999
#fit_hi_p <- permutest(fit_hi, exact=F) # still non-sig

# get confidence intervals
confint(fit_hi)

# forest plot (effects)
forest(fit_hi, slab = paste(meta_data_hi$Author, meta_data_hi$Year, sep = ", "),
       ilab = cbind(meta_data_hi$SSD.N, meta_data_hi$ASD.N), ilab.xpos = c(-9.5, -6), cex = 0.75)

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
eff_rmet <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
                  sd2i=O1.ASD.SD, data = meta_data_rmet, vtype = "LS", append = FALSE, var.names=c("eff_size","var"))

fit_rmet <- rma(yi=eff_size, vi=var, data = eff_rmet)

# permutation testing
#setseed=999
#fit_rmet_p <- permutest(fit_rmet, exact=F) # still non-sig

# get confidence intervals
confint(fit_rmet)

# forest plot (effects)
forest(fit_rmet, slab = paste(meta_data_rmet$Author, meta_data_rmet$Year, sep = ", "),
       ilab = cbind(meta_data_rmet$SSD.N, meta_data_rmet$ASD.N), ilab.xpos = c(-9.5, -6), cex = 0.75)

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


