
library(metafor)
library(plyr)
library(metaviz)

# read in data for soc cog dimensions
low_meta_data <- read.csv("/projects/loliver/Systematic_Review/Lower_Data_2020-02-24.csv")
hi_meta_data <- read.csv("/projects/loliver/Systematic_Review/Higher_Data_2020-02-24.csv")
rmet_meta_data <- read.csv("/projects/loliver/Systematic_Review/RMET_Data_2020-02-24.csv")

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

# fit model
# Restricted maximum-likelihood (REML) estimation is used by default when estimating τ2 
# the REML estimator is approximately unbiased and quite efficient; see Viechtbauer 2005; Veroniki et al., 2015
# DL (adequate and most commonly used) and PM (most recommended by Veroniki et al.) estimators produce extremely similar results 
fit_low <- rma(yi=eff_size, vi=var, data=eff_low, method="REML")

#setseed=999
#permutation testing (non-normality of observed effects)
#fit_low_p <- permutest(fit_low, exact=F)

# get confidence intervals
confint(fit_low)

# forest plot (effects) # 600 width
forest(fit_low, slab = paste(eff_low$Author, eff_low$Year, sep = ", "), xlab="Hedges' g",
         mlab="Summary (Random-Effects)", order="obs")
text(-5, 18, "Author(s) and Year",  pos=4)
text(6, 18, "Hedges' g [95% CI]", pos=2)

# if want to add Ns to forest plot
#study_table=eff_low[,c("Author","Task.SSD.N","Task.ASD.N")],table_headers=c("Author","N SSD", "N ASD"))


# outlier and influential cases check
inf_low <- influence(fit_low)
plot(inf_low, plotdfb = TRUE)

# look at externally standardized residuals in particular
rstudent(fit_low)

# run leave one out analysis
fit_low_leave1 <- leave1out(fit_low)
#write.csv(fit_low_leave1, file="/projects/loliver/Systematic_Review/Paper/Tables/fit_low_leave1.csv",row.names=F)

# generate GOSH plot - cool visualization for outliers, but takes awhile
#fit_low_gosh <- gosh(fit_low)
#plot(fit_low_gosh, out=1)

# rerun meta-analysis without outlier (study 1)
eff_low_out <- eff_low[-1,]
fit_low_out <- rma(yi=eff_size, vi=var, data=eff_low_out, method="REML")

#setseed=999
#permutation testing (non-normality of observed effects)
#fit_low_out_p <- permutest(fit_low_out, exact=T)

# get confidence intervals
confint(fit_low_out)

# forest plot
forest(fit_low_out, slab = paste(eff_low_out$Author, eff_low_out$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-3, 17.5, "Author(s) and Year",  pos=1)
text(4, 17.5, "Hedges' g [95% CI]", pos=1)

# Q-Q plot to check normality
qqnorm(fit_low_out,label="out",main = "Random-Effects Model")

# publication bias tests 
# funnel plots # 500 width
funnel(fit_low, main = "Random-Effects Model")
funnel(fit_low_out, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry
# One may be able to detect such asymmetry by testing whether the observed outcomes (or residuals from a model with moderators) 
# are related to their corresponding sampling variances, standard errors (the default, "sei"), or more simply, sample sizes
regtest(fit_low, model="rma", predictor = "sei")
regtest(fit_low_out, model="rma", predictor = "sei")

# Ns per meta
colSums(eff_low[,c("Task.SSD.N","Task.ASD.N")])
colSums(eff_low_out[,c("Task.SSD.N","Task.ASD.N")])


# moderator analyses
# with pub year as moderator (check)
fit_low_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_low)
fit_low_yr_out <- rma(yi=eff_size, vi=var, mods=Year, data=eff_low_out)

# with pub year and age as moderators (use)
fit_low_yrage <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_low)
fit_low_yrage_out <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_low_out)

# check age as categorical moderator (sig diff yes or no)
fit_low_agecat_out <- rma(yi=eff_size, vi=var, mods=~factor(Case.Age.Sig.Diff), data=eff_low_out[eff_low_out$Case.Age.Sig.Diff!="Skip",])

# check age in those with sig diff only (and those without sig diff)
rma(yi=eff_size, vi=var, mods=age_diff, data=eff_low_out[eff_low_out$Case.Age.Sig.Diff=="Yes",])
rma(yi=eff_size, vi=var, mods=age_diff, data=eff_low_out[eff_low_out$Case.Age.Sig.Diff=="No",])

# exploratory QA Score as moderator 
fit_low_qa <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_low)
fit_low_qa_out <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_low_out)

# exploratory - prop on antipsychotics - Bolte and Poustka (Study 1) not included regardless
fit_low_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_low_out)

# exploratory - prop male (don't use)
fit_low_sex <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_low)
fit_low_sex_out <- rma(yi=eff_size, vi=var, mods=Diff.Male, data=eff_low_out)


# sensitivity analyses - data availability
fit_low_sens_avail <- rma(yi=eff_size, vi=var, data=eff_low[eff_low$Data.Availability=="Yes",])
fit_low_sens_avail_out <- rma(yi=eff_size, vi=var, data=eff_low_out[eff_low_out$Data.Availability=="Yes",])

regtest(fit_low_sens_avail_out, model="rma", predictor = "sei")
colSums(eff_low_out[eff_low_out$Data.Availability=="Yes",c("Task.SSD.N","Task.ASD.N")])

# excluding early psychosis (Pepper) SSD groups (no SPD for lower-level)
fit_low_sens_SSD <-  rma(yi=eff_size, vi=var, data=eff_low[eff_low$Author!="Pepper et al.",])
fit_low_sens_SSD_out <-  rma(yi=eff_size, vi=var, data=eff_low_out[eff_low_out$Author!="Pepper et al.",])

regtest(fit_low_sens_SSD_out, model="rma", predictor = "sei")
colSums(eff_low_out[eff_low_out$Author!="Pepper et al.",c("Task.SSD.N","Task.ASD.N")])

# excluding adolescent (Bolte & Poustka, Waris; Bolte and Poustka is the outlier) studies - no children only included
fit_low_sens_young_out <-  rma(yi=eff_size, vi=var, data=eff_low_out[eff_low_out$Author!="Waris et al.",])

regtest(fit_low_sens_young_out, model="rma", predictor = "sei")
colSums(eff_low_out[eff_low_out$Author!="Waris et al.",c("Task.SSD.N","Task.ASD.N")])


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
text(-10.5, 19.75, "Author(s) and Year",  pos=1)
text(6.75, 19.75, "Hedges' g [95% CI]", pos=1)

# outlier and influential cases check
inf_hi <- influence(fit_hi)
plot(inf_hi, plotdfb = TRUE)

# look at externally standardized residuals in particular
rstudent(fit_hi)

# run leave one out analysis
fit_hi_leave1 <- leave1out(fit_hi)
#write.csv(fit_hi_leave1, file="/projects/loliver/Systematic_Review/Paper/Tables/fit_hi_leave1.csv",row.names=F)

# generate GOSH plot - cool visualization for outliers, but takes forever 
#fit_hi_gosh <- gosh(fit_hi)
#plot(fit_hi_gosh, out=10)

# rerun meta-analysis without outlier (study 10, Waris)
eff_hi_out <- eff_hi[-10,]
fit_hi_out <- rma(yi=eff_size, vi=var, data=eff_hi_out, method="REML")

# get confidence intervals
confint(fit_hi_out)

# forest plot (effects)
forest(fit_hi_out, slab = paste(eff_hi_out$Author, eff_hi_out$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-3.75, 18.5, "Author(s) and Year",  pos=1)
text(3.75, 18.5, "Hedges' g [95% CI]", pos=1)

# Q-Q plot to check normality
qqnorm(fit_hi_out,label="out",main = "Random-Effects Model")

# publication bias tests
# funnel plots
funnel(fit_hi, main = "Random-Effects Model")
funnel(fit_hi_out, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry
regtest(fit_hi, model="rma", predictor = "sei")
regtest(fit_hi_out, model="rma", predictor = "sei")

# Ns per meta
colSums(eff_hi[,c("Task.SSD.N","Task.ASD.N")])
colSums(eff_hi_out[,c("Task.SSD.N","Task.ASD.N")])


# moderator analyses
# with pub year as moderator (check)
fit_hi_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_hi)
fit_hi_yr_out <- rma(yi=eff_size, vi=var, mods=Year, data=eff_hi_out)

# with pub year and age as moderators (use)
fit_hi_yrage <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_hi)
fit_hi_yrage_out <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_hi_out)

# check age as categorical moderator (sig diff yes or no)
fit_hi_agecat_out <- rma(yi=eff_size, vi=var, mods=~factor(Case.Age.Sig.Diff), data=eff_hi_out[eff_hi_out$Case.Age.Sig.Diff!="Skip",])

# check age in those with sig diff only (and those without sig diff)
rma(yi=eff_size, vi=var, mods=age_diff, data=eff_hi_out[eff_hi_out$Case.Age.Sig.Diff=="Yes",])
rma(yi=eff_size, vi=var, mods=age_diff, data=eff_hi_out[eff_hi_out$Case.Age.Sig.Diff=="No",])

# with QA Score as moderator 
fit_hi_qa <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_hi)
fit_hi_qa_out <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_hi_out)

# stim type (verb and vis only); Waris not included either way
fit_hi_stim <- rma(yi=eff_size, vi=var, mods=~factor(Stim.Type), data=eff_hi_out[eff_hi_out$Stim.Type!="both",])

# exploratory - prop on antipsychotics - k = 8; Waris not included either way and doesn't change with update
fit_hi_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_hi)

# check to see if either group driving this
fit_hi_med_SSD <- rma(yi=eff_size, vi=var, mods=SSD.Prop.on.Antipsych, data=eff_hi)
fit_hi_med_ASD <- rma(yi=eff_size, vi=var, mods=ASD.Prop.on.Antipsych, data=eff_hi)

forest(fit_hi_med,order="obs")

# Scatterplot showing effect sizes of the individual studies plotted against diff in prop on antipsychotics 
# The radius of the points is drawn proportional to the inverse of the standard errors 
# i.e., larger/more precise studies are shown as larger points
# Hence, the area of the points is drawn proportional to the inverse sampling variances. 

size <- 1/sqrt(eff_hi$var)
size <- 3*(size/max(size))

plot(eff_hi$Diff.on.Antipsych,eff_hi$eff_size,pch=19,cex=size,
    xlab="Difference in Proportion of Participants on Antipsychotics (SSDs-ASD)",
    ylab="Hedges' g",xlim=c(0.4,0.9),ylim=c(-1,1))

abline(h=0,lty="dotted")
abline(lm(eff_hi$eff_size~eff_hi$Diff.on.Antipsych))


# sensitivity analyses - data availability
fit_hi_sens_avail <- rma(yi=eff_size, vi=var, data=eff_hi[eff_hi$Data.Availability=="Yes",])
fit_hi_sens_avail_out <- rma(yi=eff_size, vi=var, data=eff_hi_out[eff_hi_out$Data.Availability=="Yes",])

regtest(fit_hi_sens_avail_out, model="rma", predictor = "sei")
colSums(eff_hi_out[eff_hi_out$Data.Availability=="Yes",c("Task.SSD.N","Task.ASD.N")])

# excluding SPD (Booules-Katri) and early psychosis (Pepper) SSD groups
fit_hi_sens_SSD <- rma(yi=eff_size, vi=var, data=eff_hi[eff_hi$Author!="Booules-Katri et al."&eff_hi$Author!="Pepper et al.",])
fit_hi_sens_SSD_out <- rma(yi=eff_size, vi=var, data=eff_hi_out[eff_hi_out$Author!="Booules-Katri et al."&eff_hi_out$Author!="Pepper et al.",])

regtest(fit_hi_sens_SSD_out, model="rma", predictor = "sei")
colSums(eff_hi_out[eff_hi_out$Author!="Booules-Katri et al."&eff_hi_out$Author!="Pepper et al.",c("Task.SSD.N","Task.ASD.N")])

# excluding child (Pilowsky) and adolescent (Tin, Waris) studies - Waris already excluded with outlier removed
fit_hi_sens_young_out <-  rma(yi=eff_size, vi=var, data=eff_hi_out[eff_hi_out$Author!="Pilowsky et al."&eff_hi_out$Author!="Tin et al.",])

regtest(fit_hi_sens_young_out, model="rma", predictor = "sei")
colSums(eff_hi_out[eff_hi_out$Author!="Pilowsky et al."&eff_hi_out$Author!="Tin et al.",c("Task.SSD.N","Task.ASD.N")])


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
#fit_rmet_p <- permutest(fit_rmet, exact=F)

# get confidence intervals
confint(fit_rmet)

# forest plot (effects)
forest(fit_rmet, slab = paste(eff_rmet$Author, eff_rmet$Year, sep = ", "), xlab="Hedges' g",
       mlab="Summary (Random-Effects)", order="obs")
text(-3.5, 14.5, "Author(s) and Year",  pos=1)
text(4.3, 14.5, "Hedges' g [95% CI]", pos=1)

# outlier and influential cases check
inf_rmet <- influence(fit_rmet)
plot(inf_rmet, plotdfb = TRUE)

# look at externally standardized residuals in particular
rstudent(fit_rmet)

# run leave one out analysis
fit_rmet_leave1 <- leave1out(fit_rmet)
#write.csv(fit_rmet_leave1, file="/projects/loliver/Systematic_Review/Paper/Tables/fit_rmet_leave1.csv",row.names=F)

# generate GOSH plot
#fit_rmet_gosh <- gosh(fit_rmet)
#plot(fit_rmet_gosh, col="blue")

# Q-Q plot to check normality
qqnorm(fit_rmet,label="out",main = "Random-Effects Model")

# publication bias tests
# funnel plots
funnel(fit_rmet, main = "Random-Effects Model")

# Egger's regression test for funnel plot assymetry
regtest(fit_rmet, model="rma", predictor = "sei")

# Ns per meta
colSums(eff_rmet[,c("Task.SSD.N","Task.ASD.N")])


# moderator analyses
# with pub year as moderator (check)
fit_rmet_yr <- rma(yi=eff_size, vi=var, mods=Year, data=eff_rmet)

# with pub year and age as moderators (use) - p=.069 for Year (but p=.098 overall)
fit_rmet_yrage <- rma(yi=eff_size, vi=var, mods=cbind(Year,age_diff), data=eff_rmet)

# check age as categorical moderator (sig diff yes or no)
fit_rmet_agecat <- rma(yi=eff_size, vi=var, mods=~factor(Case.Age.Sig.Diff), data=eff_rmet[eff_rmet$Case.Age.Sig.Diff!="Skip",])

# check age in those with sig diff only (and those without sig diff)
rma(yi=eff_size, vi=var, mods=age_diff, data=eff_rmet[eff_rmet$Case.Age.Sig.Diff=="Yes",])
rma(yi=eff_size, vi=var, mods=age_diff, data=eff_rmet[eff_rmet$Case.Age.Sig.Diff=="No",])

# with QA Score as moderator 
fit_rmet_qa <- rma(yi=eff_size, vi=var, mods=QA.Score, data=eff_rmet)

# exploratory - prop on antipsychotics
fit_rmet_med <- rma(yi=eff_size, vi=var, mods=Diff.on.Antipsych, data=eff_rmet)


# sensitivity analyses - data availability
fit_rmet_sens_avail <- rma(yi=eff_size, vi=var, data=eff_rmet[eff_rmet$Data.Availability=="Yes",])

regtest(fit_rmet_sens_avail, model="rma", predictor = "sei")
colSums(eff_rmet[eff_rmet$Data.Availability=="Yes",c("Task.SSD.N","Task.ASD.N")])

# forest plot (effects)
forest(fit_rmet_sens_avail, slab = paste(eff_rmet[eff_rmet$Data.Availability=="Yes","Author"], 
    eff_rmet[eff_rmet$Data.Availability=="Yes","Year"], sep = ", "), xlab="Hedges' g",mlab="Summary (Random-Effects)", order="obs")
text(-2.75, 13.5, "Author(s) and Year",  pos=1)
text(4.1, 13.5, "Hedges' g [95% CI]", pos=1)

# excluding SPD (Booules-Katri) and early psychosis (Pepper) SSD groups
fit_rmet_sens_SSD <- rma(yi=eff_size, vi=var, data=eff_rmet[eff_rmet$Author!="Booules-Katri et al."&eff_rmet$Author!="Pepper et al.",])

regtest(fit_rmet_sens_SSD, model="rma", predictor = "sei")
colSums(eff_rmet[eff_rmet$Author!="Booules-Katri et al."&eff_rmet$Author!="Pepper et al.",c("Task.SSD.N","Task.ASD.N")])

# no child and adolescent only studies included


# Overall summary stats for qualitative synthesis

# read in all data extracted
all_meta_data <- read.csv("/projects/loliver/Systematic_Review/All_Data_2020-02-24.csv")

# Ns per group (total)
colSums(all_meta_data[,c("SSD.N","ASD.N")])

# excluding those not in meta
colSums(all_meta_data[c(1:8,10:12,14:16,18:34),c("SSD.N","ASD.N")])

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

