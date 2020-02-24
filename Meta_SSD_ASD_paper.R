
library(metafor)
library(plyr)
library(metaviz)

# read in data for soc cog dimensions
low_meta_data <- read.csv("/projects/loliver/Systematic_Review/Lower_Data_2020-01-09.csv")
hi_meta_data <- read.csv("/projects/loliver/Systematic_Review/Higher_Data_2020-01-09.csv")
rmet_meta_data <- read.csv("/projects/loliver/Systematic_Review/RMET_Data_2020-01-09.csv")

# calculate individual effect sizes
# specifically, calculate standardized mean effect sizes (Hedge's g) from the means, SDs, and Ns for the SSD and ASD groups for each paper
# vtype="LS" (the default) uses the usual large-sample approximation to compute the sampling variances
# vtype="UB" provides unbiased estimates of the sampling variances

# low (no RMET)
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

#setseed=999
#permutation testing (non-normality of observed effects)
#fit_low_p <- permutest(fit_low, exact=T) # still sig

# get confidence intervals
confint(fit_low)

# forest plot (effects) # 600 width
forest(fit_low, slab = paste(eff_low$Author, eff_low$Year, sep = ", "), xlab="Hedges' g",
         mlab="Summary (Random-Effects)", order="obs")
text(-4, 15, "Author(s) and Year",  pos=4)
text( 5.5, 15, "Hedges' g [95% CI]", pos=2)

# width 1250, height 500
viz_forest(x=fit_low, study_labels=eff_low$Author, summary_label = "Summary (Random-Effects)", xlab = "Hedges' g", col="Black",
           text_size=7, variant="classic", annotate_CI=T)

viz_forest(x=fit_low, study_labels=eff_low$Author, summary_label = "Summary (Random-Effects)", xlab = "Hedges' g", col="Greys",
           text_size=5, variant="rain")

# if want to add Ns to forest plot
#study_table=eff_low[,c("Author","Task.SSD.N","Task.ASD.N")],table_headers=c("Author","N SSD", "N ASD"))


# outlier and influential cases check
inf_low <- influence(fit_low)
plot(inf_low, plotdfb = TRUE)

# look at externally standardized residuals in particular
rstudent(fit_low)

# can also run leave one out analysis - stil sig with Bolte & Poustka out
leave1out(fit_low)

# generate GOSH plot
fit_low_gosh <- gosh(fit_low)
plot(fit_low_gosh, out=1)

# rerun meta-analysis without outlier (study 1)
eff_low_out <- eff_low[-1,]
fit_low_out <- rma(yi=eff_size, vi=var, data=eff_low_out, method="REML")

#setseed=999
#permutation testing (non-normality of observed effects)
#fit_low_out_p <- permutest(fit_low_out, exact=T) # still sig

# get confidence intervals
confint(fit_low_out)

# forest plot (effects)
forest(fit_low_out, slab = paste(eff_low_out$Author, eff_low_out$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-2.5, 15, "Author(s) and Year",  pos=1)
text(3.5, 15, "Hedges' g [95% CI]", pos=1)

# Q-Q plot to check normality
qqnorm(fit_low_out,label="out",main = "Random-Effects Model")

# publication bias tests
# funnel plots
funnel(fit_low, main = "Random-Effects Model")
funnel(fit_low_out, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry (pub bias)
# One may be able to detect such asymmetry by testing whether the observed outcomes (or residuals from a model with moderators) 
# are related to their corresponding sampling variances, standard errors (the default, "sei"), or more simply, sample sizes
regtest(fit_low, model="rma", predictor = "sei")
regtest(fit_low_out, model="rma", predictor = "sei")

# Ns per meta
colSums(eff_low[,c("Task.SSD.N","Task.ASD.N")])
colSums(eff_low_out[,c("Task.SSD.N","Task.ASD.N")])


# moderator analyses
# with pub year as moderator - sig (accounts for 70% of heterogeneity) - NS once outlier removed
fit_low_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_low)
fit_low_yr_out <- rma(yi=eff_size, vi=var, mods=Year, data=eff_low_out)

# with pub year and age as moderators - sig (accounts for 60% of heterogeneity) - both NS once outlier removed
fit_low_yrage <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_low)
fit_low_yrage_out <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_low_out)

# with QA Score and age as moderators
# did check age and QA alone and NS, as with together
fit_low_qaage <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_low)
fit_low_qaage_out <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_low_out)

# with QA Score as moderator 
fit_low_qa <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_low)
fit_low_qa_out <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_low_out)

# exploratory - prop on antipsychotics - NS - Bolte and Poustka (Study 1) not included regardless
fit_low_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_low)

# exploratory - IQ (full-scale or WRAT) k = 4 (don't use)
fit_low_IQ <- rma(yi=eff_size, vi=var, mods=IQ_diff, data=eff_low[eff_low$IQ.Type=="Full-scale"|eff_low$IQ.Type=="WRAT",])

# exploratory - prop male (don't use) - both NS
fit_low_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_low)
fit_low_sex_out <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_low_out)


# sensitivity analyses - data availability - NS but use
fit_low_avail <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_low)
fit_low_avail_out <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_low_out)

# excluding early psychosis (Pepper) SSD groups (no SPD for lower-level)
fit_low_sens_SSD <-  rma(yi=eff_size, vi=var, data=eff_low[c(1:11,13),])
fit_low_sens_SSD_out <-  rma(yi=eff_size, vi=var, data=eff_low_out[c(1:10,12),])

# excluding PDD ASD group (Waris)
fit_low_sens_ASD <-  rma(yi=eff_size, vi=var, data=eff_low[c(1:9,11:13),])
fit_low_sens_ASD_out <-  rma(yi=eff_size, vi=var, data=eff_low_out[c(1:8,10:12),])

# excluding adolescent (Bolte & Poustka, Waris) studies - no children only included - skip
fit_low_sens_young <-  rma(yi=eff_size, vi=var, data=eff_low[c(2:9,11:13),]) # same as above with outlier removed


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

# permutation testing
#setseed=999
#fit_hi_p <- permutest(fit_hi, exact=F) # still non-sig

# get confidence intervals
confint(fit_hi)

# forest plot (effects)
forest(fit_hi, slab = paste(eff_hi$Author, eff_hi$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-11, 19, "Author(s) and Year",  pos=1)
text(6.75, 19, "Hedges' g [95% CI]", pos=1)

# outlier and influential cases check
inf_hi <- influence(fit_hi)
plot(inf_hi, plotdfb = TRUE)

# look at externally standardized residuals in particular
rstudent(fit_hi)

# can also run leave one out analysis
leave1out(fit_hi)

# generate GOSH plot
fit_hi_gosh <- gosh(fit_hi)
plot(fit_hi_gosh, out=10)

# rerun meta-analysis without outlier (study 1)
eff_hi_out <- eff_hi[-10,]
fit_hi_out <- rma(yi=eff_size, vi=var, data=eff_hi_out, method="REML")

# get confidence intervals
confint(fit_hi_out)

# forest plot (effects)
forest(fit_hi_out, slab = paste(eff_hi_out$Author, eff_hi_out$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-3.75, 18, "Author(s) and Year",  pos=1)
text(3.75, 18, "Hedges' g [95% CI]", pos=1)

# Q-Q plot to check normality
qqnorm(fit_hi_out,label="out",main = "Random-Effects Model")

# publication bias tests
# funnel plots
funnel(fit_hi, main = "Random-Effects Model")
funnel(fit_hi_out, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry (pub bias)
regtest(fit_hi, model="rma", predictor = "sei") # sig
regtest(fit_hi_out, model="rma", predictor = "sei") # NS

# Ns per meta
colSums(eff_hi[,c("Task.SSD.N","Task.ASD.N")])
colSums(eff_hi_out[,c("Task.SSD.N","Task.ASD.N")])


# moderator analyses
# with pub year as moderator
fit_hi_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_hi)
fit_hi_yr_out <- rma(yi=eff_size, vi=var, mods=Year, data=eff_hi_out)

# with pub year and age as moderators
fit_hi_yrage <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_hi)
fit_hi_yrage_out <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_hi_out)

# with QA Score and age as moderators
fit_hi_qaage <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_hi)
fit_hi_qaage_out <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_hi_out)

# with QA Score as moderator 
fit_hi_qa <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_hi)
fit_hi_qa_out <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_hi_out)

# stim type (verb and vis only) - NS; Waris not included either way
fit_hi_stim <- rma(yi=eff_size, vi=var, mods=~factor(Stim.Type), data=eff_hi[eff_hi$Stim.Type!="both",])

# exploratory - prop on antipsychotics - sig, k = 8; Waris not included either way
fit_hi_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_hi)

# check to see if either group driving this
fit_hi_med_SSD <- rma(yi=eff_size, vi=var, mods=SSD.Prop.on.Antipsych, data=eff_hi)
fit_hi_med_ASD <- rma(yi=eff_size, vi=var, mods=ASD.Prop.on.Antipsych, data=eff_hi)

forest(fit_hi_med,order="obs")

# calculate predicted effect sizes for .4-.9 antipsychotic diff 
preds_med <- predict(fit_hi_med,newmods=seq(from=0.4,to=.9,by=.05))

# area of the points will be proportional to inverse variances (from metafor example)
size <- 1 / sqrt(eff_hi$var)
size <- 3*(size / max(size))

plot(eff_hi$Diff.on.Antipsych,eff_hi$eff_size,pch=19,cex=size,
    xlab="Difference in Proportion of Participants on Antipsychotics (SSD-ASD)",
    ylab="Hedges' g",xlim=c(0.4,0.9),ylim=c(-1,1))
 
#text(eff_hi$Diff.on.Antipsych,eff_hi$eff_size,eff_hi[!is.na(eff_hi$Diff.on.Antipsych),"Author"],adj=c(0,2))
abline(h=0,lty="dotted")
abline(lm(eff_hi$eff_size~eff_hi$Diff.on.Antipsych))

# add predicted values and corresponding CI bounds - unnecessary I think
# lines(seq(from=0.4,to=.9,by=.05), preds_med$pred)
# lines(seq(from=0.4,to=.9,by=.05), preds_med$ci.lb,lty="dashed")
# lines(seq(from=0.4,to=.9,by=.05), preds_med$ci.ub,lty="dashed")


# exploratory - IQ (full-scale or WRAT); Waris not included either way (don't use)
fit_hi_IQ <- rma(yi=eff_size, vi=var, mods=IQ_diff, data=eff_hi[eff_hi$IQ.Type=="Full-scale"|eff_hi$IQ.Type=="WRAT",])

# exploratory - prop male - sig without Waris (don't use)
fit_hi_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_hi)
fit_hi_sex_out <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_hi_out)


# sensitivity analyses - data availability - all NS
fit_hi_avail <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_hi)
fit_hi_avail_out <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_hi_out)

# excluding SPD (Booules-Katri) and early psychosis (Pepper) SSD groups
fit_hi_sens_SSD <- rma(yi=eff_size, vi=var, data=eff_hi[c(1:10,12,14:16),])
fit_hi_sens_SSD_out <- rma(yi=eff_size, vi=var, data=eff_hi_out[c(1:9,11,13:15),])

# excluding PDD ASD group (Waris) - already done as this is the outlier - skip
fit_hi_sens_ASD <-  rma(yi=eff_size, vi=var, data=eff_hi[-10,])

# excluding child (Pilowsky) and adolescent (Tin, Waris) studies - Waris already excluded with outlier removed
fit_hi_sens_young <-  rma(yi=eff_size, vi=var, data=eff_hi[c(1:6,8,11:16),])


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

# permutation testing
#setseed=999
#fit_rmet_p <- permutest(fit_rmet, exact=F) # still non-sig

# get confidence intervals
confint(fit_rmet)

# forest plot (effects)
forest(fit_rmet, slab = paste(eff_rmet$Author, eff_rmet$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-3.75, 13.75, "Author(s) and Year",  pos=1)
text(4.3, 13.75, "Hedges' g [95% CI]", pos=1)

# outlier and influential cases check
inf_rmet <- influence(fit_rmet)
plot(inf_rmet, plotdfb = TRUE)

# look at externally standardized residuals in particular
rstudent(fit_rmet)

# can also run leave one out analysis - all NS
leave1out(fit_rmet)

# generate GOSH plot
fit_rmet_gosh <- gosh(fit_rmet)
plot(fit_rmet_gosh, col="blue")

# Q-Q plot to check normality
qqnorm(fit_rmet,label="out",main = "Random-Effects Model")

# publication bias tests
# funnel plots
funnel(fit_rmet, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry (pub bias)
# One may be able to detect such asymmetry by testing whether the observed outcomes (or residuals from a model with moderators) 
# are related to their corresponding sampling variances, standard errors (the default, "sei"), or more simply, sample sizes
regtest(fit_rmet, model="rma", predictor = "sei")

# Ns per meta
colSums(eff_rmet[,c("Task.SSD.N","Task.ASD.N")])


# moderator analyses
# with pub year as moderator
fit_rmet_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_rmet)

# with pub year and age as moderators
fit_rmet_yrage <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_rmet)

# with QA Score and age as moderators
fit_rmet_qaage <- rma(yi=eff_size, vi=var, mods=cbind(QA.Score,age_diff), data=eff_rmet)

# with QA Score as moderator 
fit_rmet_qa <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_rmet)

# exploratory - prop on antipsychotics
fit_rmet_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_rmet)

# exploratory - IQ (full-scale or WRAT)
fit_rmet_IQ <- rma(yi=eff_size, vi=var, mods=IQ_diff, data=eff_rmet[eff_rmet$IQ.Type=="Full-scale"|eff_rmet$IQ.Type=="WRAT",])

# exploratory - prop male
fit_rmet_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_rmet)


# sensitivity analyses - data availability - NS but use
fit_rmet_avail <- rma(yi=eff_size, vi=var, mods=~factor(Data.Availability), data=eff_rmet)

# excluding SPD (Booules-Katri) and early psychosis (Pepper) SSD groups
fit_rmet_sens_SSD <-  rma(yi=eff_size, vi=var, data=eff_rmet[c(1:7,9,11),])

# no PDD ASD group (Waris) included

# no child and adolescent only studies included


# Overall summary stats for qualitative synthesis

# read in all data extracted
all_meta_data <- read.csv("/projects/loliver/Systematic_Review/All_Data_2020-02-04.csv")

# Ns per group
colSums(all_meta_data[,c("SSD.N","ASD.N")])

# without those not in meta
colSums(all_meta_data[c(1:7,9:12,14:16,18:30),c("SSD.N","ASD.N")])

# average prop male
colMeans(all_meta_data[,c("SSD.Prop.Male","ASD.Prop.Male")],na.rm=T)

# studies with only males
all_meta_data[all_meta_data$SSD.Prop.Male==1&all_meta_data$ASD.Prop.Male==1,c("Author","Year","SSD.Prop.Male","ASD.Prop.Male")]

# age per group

## N: vector of sizes
## M: vector of means
## S: vector of standard deviations

grand.mean <- function(M, N) {weighted.mean(M, N)}
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}

grand.mean(all_meta_data[!is.na(all_meta_data$SSD.Age.Mean),"SSD.Age.Mean"],all_meta_data[!is.na(all_meta_data$SSD.Age.Mean),"SSD.N"])
grand.sd(all_meta_data[!is.na(all_meta_data$SSD.Age.Mean),"SSD.Age.SD"],all_meta_data[!is.na(all_meta_data$SSD.Age.Mean),"SSD.Age.Mean"],all_meta_data[!is.na(all_meta_data$SSD.Age.Mean),"SSD.N"])

grand.mean(all_meta_data[!is.na(all_meta_data$ASD.Age.Mean),"ASD.Age.Mean"],all_meta_data[!is.na(all_meta_data$ASD.Age.Mean),"ASD.N"])
grand.sd(all_meta_data[!is.na(all_meta_data$ASD.Age.Mean),"ASD.Age.SD"],all_meta_data[!is.na(all_meta_data$ASD.Age.Mean),"ASD.Age.Mean"],all_meta_data[!is.na(all_meta_data$ASD.Age.Mean),"ASD.N"])


# IQ type
count(all_meta_data$IQ.Type)

# Mean full-scale IQ
grand.mean(all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$SSD.IQ.Mean),"SSD.IQ.Mean"],all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$SSD.IQ.Mean),"SSD.N"])
grand.sd(all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$SSD.IQ.Mean),"SSD.IQ.SD"],all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$SSD.IQ.Mean),"SSD.IQ.Mean"],
         all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$SSD.IQ.Mean),"SSD.N"])

grand.mean(all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$ASD.IQ.Mean),"ASD.IQ.Mean"],all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$ASD.IQ.Mean),"ASD.N"])
grand.sd(all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$ASD.IQ.Mean),"ASD.IQ.SD"],all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$ASD.IQ.Mean),"ASD.IQ.Mean"],
         all_meta_data[all_meta_data$IQ.Type=="Full-scale"&!is.na(all_meta_data$ASD.IQ.Mean),"ASD.N"])

# sample type
all_meta_data[,c("Sample","SSD.Sample","ASD.Sample")]

# CPZ-eq
all_meta_data[,c("SSD.Prop.on.Antipsych","SSD.Mean.CPZ.eq..mg.","ASD.Prop.on.Antipsych","ASD.Mean.CPZ.eq..mg.")]

