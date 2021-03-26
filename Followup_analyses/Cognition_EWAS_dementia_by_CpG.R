##################################################################################################################
# Name: Cognition_EWAS_dementia_by_CpG

# Purpose: Run linear regression analyses between DNA methylation and dementia status (to make sure associations
# are not just a result of cognitive decline in preclinical dementia)
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20210311 by Ida Karlsson (IK) 
# Updated:	

##################################################################################################################
##################################################################################################################

.libPaths("Z:/Programs/R/Packages")
library(multiwayvcov)
library(lmtest)

# Read in meth data
adata <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20191216.txt", header=T)
adata1 <- adata[c("TWINNR", "cg04549090", "cg18064256")]
adata2 <- adata[,2:9]

adata3 <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tmemod_20191216.txt", header=T)
adata4 <- adata3[c("TWINNR", "cg09988380", "cg25651129", "cg08011941")]

adata5 <- merge(adata2, adata1, by="TWINNR")
adata6 <- merge(adata5, adata4, by="TWINNR")

# Read in dementia data
dem <- read.table("P:/Dementia_IK/Various data across and outside studies/Master_dementia_file/Data/DerivedData/master_dementia_20191107.txt", header=T, sep="\t")
dem <- dem[c(1,15)]

# Merge with cog data
dem_adata <- merge(adata6, dem, by.x="TWINNR", by.y="twinnr")

# Linear models (done one by one)
lr_CpG_adj <- lm(cg04549090 ~ gosh_dem + AGE + SEX + SMOKE + CHIP, data=dem_adata)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(dem_adata$PAIRID)) 
coefficients1 <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients1
#               Estimate  Std. Error t value Pr(>|t|)
#gosh_dem     8.4648e-02  1.2210e-01  0.6933   0.4884

lr_CpG_adj <- lm(cg18064256 ~ gosh_dem + AGE + SEX + SMOKE + CHIP, data=dem_adata)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(dem_adata$PAIRID)) 
coefficients2 <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients2
#               Estimate  Std. Error t value Pr(>|t|)
#gosh_dem     0.0180788  0.1282018  0.1410 0.887909 

lr_CpG_adj <- lm(cg09988380 ~ gosh_dem + AGE + SEX + SMOKE + CHIP, data=dem_adata)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(dem_adata$PAIRID)) 
coefficients3 <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients3
#               Estimate  Std. Error t value Pr(>|t|)
#gosh_dem    -0.0059276  0.1226388 -0.0483 0.961468

lr_CpG_adj <- lm(cg25651129 ~ gosh_dem + AGE + SEX + SMOKE + CHIP, data=dem_adata)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(dem_adata$PAIRID)) 
coefficients4 <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients4
#               Estimate  Std. Error t value Pr(>|t|)
#gosh_dem     0.0502730  0.1453560  0.3459   0.7296

lr_CpG_adj <- lm(cg08011941 ~ gosh_dem + AGE + SEX + SMOKE + CHIP, data=dem_adata)
lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(dem_adata$PAIRID)) 
coefficients5 <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
coefficients5
#               Estimate  Std. Error t value Pr(>|t|)
#gosh_dem     0.0151294  0.1348719  0.1122  0.91073 

### --> No evidence of difference in DNA methylation by dementia status

##################################################################################################################
##################################################################################################################
##################################################################################################################
