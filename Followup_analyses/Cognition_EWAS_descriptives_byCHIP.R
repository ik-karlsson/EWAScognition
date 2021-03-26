##################################################################################################################
# Name: Cognition_EWAS_descriptives_byCHIP

# Purpose: Descriptive analyses, for the full sample and by chip
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20210311 by Ida Karlsson (IK) 
# Updated:	

##################################################################################################################
##################################################################################################################
 
.libPaths("Z:/Programs/R/Packages")
library(plyr)
library(psych)
library(gmodels)

# Read in the cognitive data
cog <- read.table("P:/Dementia_IK/Cognition/Cognitive_data/Cognition_data_long_20191216.txt", header=T, sep="\t")
names(cog)[1] <- "TWINNR"
names(cog)[7] <- "cog_age"

# Read in pheno data
pdata <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20191216.txt", header=T)
pdata <- pdata[2:9]

# Read in meth data ("raw M-values" extracted)
mdata <- read.table("P:/Dementia_IK/Cognition/EWAS/Data/Original_data/Extract_CpG_Mvalues_20210324.txt", header=T)

# Merge pheno data and meth data
pmdata <- merge(pdata, mdata, by="TWINNR")

# Merge with cog data
adata <- merge(pmdata, cog, by="TWINNR")

# Count number of samples
cog_adata <- ddply(cog_adata,.(TWINNR),transform, N_COG = NROW(piece))

### First and last measure of cognition to look at follow-up
# First measure
# NOTE! Keeping cognitive measures
cog_first <- cog_adata[order(cog_adata$TWINNR, cog_adata$cog_age),]
cog_first <- cog_first[!duplicated(cog_first$TWINNR),]

cog_first <- cog_first[c(1,7,2:5,8,9)]
colnames(cog_first)[2] <- ("FIRSTAGE_COG")

# Last
# NOTE! Also keeping number of measures from last observation
cog_last <- cog_adata[order(cog_adata$TWINNR, -cog_adata$cog_age),]
cog_last <- cog_last[!duplicated(cog_last$TWINNR),]

cog_last <- cog_last[c("TWINNR","cog_age", "N_COG")]
colnames(cog_last)[2] <- ("LASTAGE_COG")

# Merge together and calculate follow-up time 
cog_fup <- merge(cog_first, cog_last, by="TWINNR")
cog_fup$COG_FUP <- cog_fup$LASTAGE_COG - cog_fup$FIRSTAGE_COG

# Merge with adata 
adata8 <- merge(pmdata, cog_fup, by="TWINNR")

# Smoking as binary
adata8$SMOKE.b <- ifelse((adata8$SMOKE == 1), 0,
                            ifelse((adata8$SMOKE == 2 ), 0,
                                   ifelse((adata8$SMOKE == 3), 1, NA)))  

# Individuals with 3+ cognitive measurse
adata8$cog3plus <-ifelse((adata8$N_COG < 3), 0,
                         ifelse((adata8$N_COG >= 3 ), 1, NA))

# Order the columns (as in the outupt table, for convenience)
adata9 <- adata8[c(1,8,25,5,24,6,13,10,9,11,12,14,23,22,16,15,17,19,20,18)]

##################################################################################################################
### Create the table

#### Output file 
desc_results <- matrix(NA, nrow = 20, ncol = 5)
desc_results[2,1] <- "Total_N"
desc_results[3:20,1] <- colnames(adata9)[3:20]
desc_results[1,1:5] <- c("", "Total sample", "450K", "EPIC", "P-value")

# Start with total N and N with 3+ measures
N <- table(adata9$CHIP)
ntot <- N[1] + N[2]
n450k <- N[1] 
nepic <- N[2]
desc_results[2,2] <- ntot
desc_results[2,3] <- n450k
desc_results[2,4] <- nepic

N3plus <- table(adata9$CHIP,adata9$cog3plus)
desc_results[3,2] <- N3plus[1,2] + N3plus[2,2]
desc_results[3,3] <- N3plus[1,2]
desc_results[3,4] <- N3plus[2,2]


# Then the (only) 2 categorical variables: sex and smoking
### Sex (women)
x2 <- chisq.test(adata9$SEX, adata9$CHIP)
etot <- x2$observed[2,1]+x2$observed[2,2]
e450k <- x2$observed[2,1]
eepic <- x2$observed[2,2]
# N and %
desc_results[4,2] <- paste(etot, " (", formatC(100*(etot/ntot), digits = 2, format = "f"), ")" )
desc_results[4,3] <- paste(e450k, " (", formatC(100*(e450k/n450k), digits = 2, format = "f"), ")" )
desc_results[4,4] <- paste(eepic, " (", formatC(100*(eepic/nepic), digits = 2, format = "f"), ")" )
# P-value
desc_results[4,5] <- formatC(x2$p.value, digits = 2, format = "f")

# Smoking
x2 <- chisq.test(adata9$SMOKE.b, adata9$CHIP)
etot <- x2$observed[2,1]+x2$observed[2,2]
e450k <- x2$observed[2,1]
eepic <- x2$observed[2,2]
# N and %
desc_results[5,2] <- paste(etot, " (", formatC(100*(etot/ntot), digits = 2, format = "f"), ")" )
desc_results[5,3] <- paste(e450k, " (", formatC(100*(e450k/n450k), digits = 2, format = "f"), ")" )
desc_results[5,4] <- paste(eepic, " (", formatC(100*(eepic/nepic), digits = 2, format = "f"), ")" )
# P-value
desc_results[5,5] <- formatC(x2$p.value, digits = 2, format = "f")

### Continous variables are done in loop
for (j in 6:20){
# Full sample
  tt1 <- describe(adata9[,j])
  desc_results[j,2] <- paste(formatC(tt1[1,3], digits = 2, format = "f"), " (", formatC(tt1[1,4], digits = 2, format = "f"), ")" )

# By chip
  tt2 <- describeBy(adata9[,j],adata9[,2],mat=TRUE,digits=2)
  desc_results[j,3] <- paste(tt2[1,5], " (", tt2[1,6], ")" )
  desc_results[j,4] <- paste(tt2[2,5], " (", tt2[2,6], ")" )

  tt3 <- t.test(adata9[,j]~adata9[,2])  
  desc_results[j,5] <- format.pval(tt3$p.value, digits=2, eps=1e-2)
  
}  

##########  Export table
write.table(desc_results, "P:/Dementia_IK/Cognition/EWAS/Output/Reviewer_comments/Cognition_EWAS_descriptives_byCHIP_20210311.txt", row.names = F, col.names = F, quote = F, sep="\t")

##################################################################################################################

### NOTE! Check for differences in effect of cg18064256 across chips, as it was significantly different.. 
### Spatial ability
adata <- as.data.frame(read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20191216.txt", header=T,stringsAsFactors = F))

adata_450k <- adata[which(adata$CHIP=='450K'),]
adata_450k <- adata_450k[c(2,3,6:8,10,11,25)]
adata_450k2 <- as.data.frame(sapply(adata_450k, as.numeric)) 

lr_CpG_adj <- lm(qI65_tspat ~ cg18064256 + AGE + SEX + SMOKE, data=adata_450k2, weights=qI65_tspatse_inv)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(adata_450k2$PAIRID)) 
coefficients <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients
#             Estimate Std. Error t value  Pr(>|t|)    
# cg18064256  -2.125509   0.455064 -4.6708 4.173e-06 ***

adata_EPIC <- adata[which(adata$CHIP=='EPIC'),]
adata_EPIC <- adata_EPIC[c(2,3,6:8,10,11,25)]
adata_EPIC2 <- as.data.frame(sapply(adata_EPIC, as.numeric)) 

lr_CpG_adj <- lm(qI65_tspat ~ cg18064256 + AGE + SEX + SMOKE, data=adata_EPIC2, weights=qI65_tspatse_inv)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(adata_EPIC2$PAIRID)) 
coefficients <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients
#             Estimate Std. Error t value  Pr(>|t|)    
# cg18064256  -1.730223   0.666860 -2.5946 0.010457 *

### Processing speed
adata <- as.data.frame(read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tsped_20191216.txt", header=T,stringsAsFactors = F))

adata_450k <- adata[which(adata$CHIP=='450K'),]
adata_450k <- adata_450k[c(2,3,6:8,10,11,33)]
adata_450k2 <- as.data.frame(sapply(adata_450k, as.numeric)) 

lr_CpG_adj <- lm(qI65_tsped ~ cg18064256 + AGE + SEX + SMOKE, data=adata_450k2, weights=qI65_tspedse_inv)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(adata_450k2$PAIRID)) 
coefficients <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients
#             Estimate Std. Error t value  Pr(>|t|)    
# cg18064256  -1.931653   0.355302 -5.4366 9.745e-08 ***

adata_EPIC <- adata[which(adata$CHIP=='EPIC'),]
adata_EPIC <- adata_EPIC[c(2,3,6:8,10,11,33)]
adata_EPIC2 <- as.data.frame(sapply(adata_EPIC, as.numeric)) 

lr_CpG_adj <- lm(qI65_tsped ~ cg18064256 + AGE + SEX + SMOKE, data=adata_EPIC2, weights=qI65_tspedse_inv)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(adata_EPIC2$PAIRID)) 
coefficients <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients
#             Estimate Std. Error t value  Pr(>|t|)    
# cg18064256  -1.354215   0.619486 -2.1860 0.030429 *

## Estimates are comparable!

##################################################################################################################
##################################################################################################################
##################################################################################################################

