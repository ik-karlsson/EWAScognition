--- 7: ##################################################################################################################
# Name: Cognition_EWAS_mQTL_fup_20200210

# Purpose: Follow-up on mQTL effects
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20200210 by Ida Karlsson (IK) 

##################################################################################################################

# Load required libraries
.libPaths("Z:/Programs/R/Packages")
library(lmtest)
library(Formula)
library(multiwayvcov)
library(lme4)
library(stringr)
library(mediation)

# Import the cognitive and methylation data
# Spatial ability
meth_spat <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20191216.txt", header=T, sep="\t")
# Working memory
meth_wmem <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tmemod_20191216.txt", header=T, sep="\t")

# Relevant SNP data
SNP_chrXf <- read.table("P:/Dementia_IK/Cognition/EWAS/Data/Original_data/t_extracted_chrX118976619Ifem_A1dosages.txt", header=F)
SNP_chrXm <- read.table("P:/Dementia_IK/Cognition/EWAS/Data/Original_data/t_extracted_chrX118976619Imale_A1dosages.txt", header=F)
SNP_rs144 <- read.table("P:/Dementia_IK/Cognition/EWAS/Data/Original_data/t_extracted_rs144382559_A1dosages.txt", header=F)
# Stack chrX females and males
SNP_chrX <- rbind(SNP_chrXf, SNP_chrXm)
# Recode to count minor allele instead of major..
SNP_chrX$SNP <- abs(2-SNP_chrX$V3)
SNP_rs144$SNP <- abs(2-SNP_rs144$V3)

# Merge cognition/methylation data with SNP data
spat_a <- merge(meth_spat, SNP_chrX, by.x="TWINNR", by.y="V2", all.x=T)
wmem_a <- merge(meth_wmem, SNP_rs144, by.x="TWINNR", by.y="V2", all.x=T)
# n=435 have both methylation and SNP data

# Can impute co-twins for some MZs
spat_a$pSNP <- na.aggregate(spat_a$SNP, by = list(spat_a$PAIRID), FUN=max, na.rm=TRUE)
# NOTE! Warning messages about no non-missing arguments - these are all for singletons with 
# non genotype data
spat_a$SNPimp <- ifelse((!is.na(spat_a$SNP)), spat_a$SNP,
                 ifelse((is.na(spat_a$SNP) & spat_a$BESTZYG==1 & spat_a$pSNP>=0), spat_a$pSNP, NA))

check <- spat_a[c("PAIRID", "BESTZYG", "SNP", "pSNP", "SNPimp")]
# Looks good

# Same for working memory
wmem_a$pSNP <- na.aggregate(wmem_a$SNP, by = list(wmem_a$PAIRID), FUN=max, na.rm=TRUE)
# NOTE! Warning messages about no non-missing arguments - these are all for singletons with 
# non genotype data
wmem_a$SNPimp <- ifelse((!is.na(wmem_a$SNP)), wmem_a$SNP,
                        ifelse((is.na(wmem_a$SNP) & wmem_a$BESTZYG==1 & wmem_a$pSNP>=0), wmem_a$pSNP, NA))
check <- wmem_a[c("PAIRID", "BESTZYG", "SNP", "pSNP", "SNPimp")]
# Looks good

# Splitting the X chromosome data by sex for analyses..
spat_aF <- spat_a[which(spat_a$SEX==2),]
spat_aM <- spat_a[which(spat_a$SEX==1),]

##################################################################################################################
########## Regression models

### Working memory
# The SNPs on methylation
lr_rs144_meth <- lm(cg08011941 ~ SNPimp + AGE + SEX + SMOKE + CHIP, data=wmem_a)
lr_rs144_meth.vcovCL<-cluster.vcov(lr_rs144_meth, as.vector(wmem_a$PAIRID)) 
coefficients <- coeftest(lr_rs144_meth, lr_rs144_meth.vcovCL)
coefficients
#              Estimate   Std. Error t value Pr(>|t|)  
# SNPimp      -0.0232223  0.3044292 -0.0763  0.93923

# The SNP on cognition
lr_rs144_cog <- lm(qI65_tmemod ~ SNPimp + AGE + SEX + SMOKE + CHIP, data=wmem_a)
lr_rs144_cog.vcovCL<-cluster.vcov(lr_rs144_cog, as.vector(wmem_a$PAIRID)) 
coefficients <- coeftest(lr_rs144_cog, lr_rs144_cog.vcovCL)
coefficients
#             Estimate    Std. Error  t value   Pr(>|t|)
# SNPimp      -1.868129   1.931375    -0.9673   0.3339

### Spatial
# The SNPs on methylation
## Women
lr_chrX_meth <- lm(cg04549090 ~ SNPimp + AGE + SEX + SMOKE + CHIP, data=spat_aF)
lr_chrX_meth.vcovCL<-cluster.vcov(lr_chrX_meth, as.vector(spat_aF$PAIRID)) 
coefficients <- coeftest(lr_chrX_meth, lr_chrX_meth.vcovCL)
coefficients
#              Estimate    Std. Error t value   Pr(>|t|)  
# SNPimp      -0.08125646  0.27068140 -0.3002   0.7642

## Men
lr_chrX_meth <- lm(cg04549090 ~ SNPimp + AGE + SEX + SMOKE + CHIP, data=spat_aM)
lr_chrX_meth.vcovCL<-cluster.vcov(lr_chrX_meth, as.vector(spat_aM$PAIRID)) 
coefficients <- coeftest(lr_chrX_meth, lr_chrX_meth.vcovCL)
coefficients
#              Estimate   Std. Error t value  Pr(>|t|)  
# SNPimp       0.1085839  0.1891631  0.5740   0.5666

# The SNP on cognition
# Women
# Convert the outcome to numeric..
spat_aF$qI65_tspat <- as.numeric(spat_aF$qI65_tspat)
lr_chrX_cog <- lm(qI65_tspat ~ SNPimp + AGE + SEX + SMOKE + CHIP, data=spat_aF)
lr_chrX_cog.vcovCL<-cluster.vcov(lr_chrX_cog, as.vector(spat_aF$PAIRID)) 
coefficients <- coeftest(lr_chrX_cog, lr_chrX_cog.vcovCL)
coefficients
#             Estimate    Std. Error t value   Pr(>|t|)
# SNPimp       26.98306   51.05575  0.5285    0.5975  

# Men
# Convert the outcome to numeric..
spat_aM$qI65_tspat <- as.numeric(spat_aM$qI65_tspat)
lr_chrX_cog <- lm(qI65_tspat ~ SNPimp + AGE + SEX + SMOKE + CHIP, data=spat_aM)
lr_chrX_cog.vcovCL<-cluster.vcov(lr_chrX_cog, as.vector(spat_aM$PAIRID)) 
coefficients <- coeftest(lr_chrX_cog, lr_chrX_cog.vcovCL)
coefficients
#             Estimate    Std. Error t value   Pr(>|t|)
# SNPimp      -16.31317   25.46928 -0.6405  0.522550

##################################################################################################################
##################################################################################################################
##################################################################################################################
