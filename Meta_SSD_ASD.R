
library(metafor)

# read in data for lower level soc cog
meta_data_low <- read.csv("/projects/loliver/Systematic_Review/Approved_Paper_Data_Lower.csv")
meta_data_hi <- read.csv("/projects/loliver/Systematic_Review/Approved_Paper_Data_Higher.csv")
meta_data_rmet <- read.csv("/projects/loliver/Systematic_Review/Approved_Paper_Data_RMET.csv")

# calculate individual effect sizes
# specifically, calculate standardized mean effect sizes (Hedge's g) from the means, SDs, and Ns for the SSD and ASD groups for each paper
# vtype="LS" (the default) uses the usual large-sample approximation to compute the sampling variances
# vtype="UB" provides unbiased estimates of the sampling variances

# low without RMET
eff_low <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
       sd2i=O1.ASD.SD, data = meta_data_low[meta_data_low$Outcome.1!="RMET",], vtype = "LS", append = FALSE, var.names=c("eff_size","var"))

# fit 
rma(yi=eff_size, vi=var, data = eff_low)


eff_hi <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
                  sd2i=O1.ASD.SD, data = meta_data_hi, vtype = "LS", append = FALSE, var.names=c("eff_size","var"))

rma(yi=eff_size, vi=var, data = eff_hi)


eff_rmet <- escalc(measure="SMD", n1i=SSD.N, n2i=ASD.N, m1i=O1.SSD.Mean, m2i=O1.ASD.Mean, sd1i=O1.SSD.SD, 
                  sd2i=O1.ASD.SD, data = meta_data_rmet, vtype = "LS", append = FALSE, var.names=c("eff_size","var"))

rma(yi=eff_size, vi=var, data = eff_rmet)

