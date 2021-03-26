##################################################################################################################
# Name: Cognition_EWAS_Plots_20190619

# Purpose: Generate plots of change in cognitive abilities for BGA poster
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20190619 by Ida Karlsson (IK) 
# Updated:	

##################################################################################################################
##### Load libraries
.libPaths("Z:/Programs/R/Packages")
library(ggplot2)
library(stringr)
library(gridExtra)
library(lmtest)
library(Formula)
library(multiwayvcov)
library(lme4)
library(dplyr)
library(gtools)

# Read in the cognitive data
cog <- read.table("P:/Dementia_IK/Cognition/Cognitive_data/Cognition_data_long_20190110.txt", header=T, sep="\t")
names(cog)[7] <- "cog_age"


# Read in the relevant methylation data
cg04767932 <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tpcom_20181217.txt", header=T, sep=" ")
cg04767932 <- cg04767932[,c(1:9, 12)]

cg18064256 <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20181217.txt", header=T, sep=" ")
cg18064256 <- cg18064256[,c("TWINNR", "cg18064256")]

meth <- merge(cg04767932, cg18064256, by="TWINNR")

# Standardize 
meth$cg04767932 <- scale(meth$cg04767932, scale = TRUE)
meth$cg18064256 <- scale(meth$cg18064256, scale = TRUE)

# Categorize into tertiles
meth$tert_cg04767932 <- factor(ntile(meth$cg04767932, 3)) 
meth$tert_cg180642562 <- factor(ntile(meth$cg18064256, 3)) 

# Merge with cognitive data
cog_adata <- merge(cog, meth, by="TWINNR")

### Modifying variables for regression models

# Smoking as binary (current smoker or non-smoker)
# Assign incident and prevalent dementia (0= no dementia, 1=prevalent, 2=incident)
cog_adata$SMOKE.b <- ifelse((cog_adata$SMOKE == 1), 0,
                            ifelse((cog_adata$SMOKE == 2 ), 0,
                                   ifelse((cog_adata$SMOKE == 3), 1, NA)))       

# Sex as categorical
cog_adata$SEX.f <- factor(cog_adata$SEX)
is.factor(cog_adata$SEX.f)

# Chip as categorical
cog_adata$CHIP.f <- factor(cog_adata$CHIP)
is.factor(cog_adata$CHIP.f)

# Create centering age
cog_adata$cage_10 <- (cog_adata$cog_age-65)/10
cog_adata$cage2_10 <- ((cog_adata$cog_age^2)-(65^2))/19

### cg04767932 on general cognition
gc_tpcom_CpG = lmer(tpcom ~ cg04767932 + cage_10 + cage2_10 + cg04767932:cage_10 +cg04767932:cage2_10 + SEX.f + SMOKE.b + CHIP.f + 
                     (cage_10 + cage2_10 |PAIRID/TWINNR), data=cog_adata, lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)), REML = F)

# Plot curve by tertiles of methylation level
cog_adata$gc_predict <-predict(gc_tpcom_CpG, na.action = na.exclude)

ggplot(data=cog_adata, aes(x=cog_age, y=tpcom, group=tert_cg04767932, colour=tert_cg04767932)) +
  stat_smooth(aes(y = gc_predict),size=1.2, method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Age", y= "General cognitive ability", title ="cg04767932 on general cognitive ability" , color = "Methylation tertile", size=2) +
  scale_color_manual(labels = c("Low", "Middle", "High"), values = c("blue4", "darkcyan", "darkorchid2")) +
  coord_cartesian(xlim=c(55, 90)) +
  theme(plot.title = element_text(size=20, face="bold"),
        axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=15),
        legend.justification = c(0, 0), legend.position = c(0, 0),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20)) 

### cg18064256 on spatial ability
gc_tspat_CpG = lmer(tspat ~ cg18064256 + cage_10 + cage2_10 + cg18064256:cage_10 +cg18064256:cage2_10 + SEX.f + SMOKE.b + CHIP.f + 
                      (cage_10 + cage2_10 |PAIRID/TWINNR), data=cog_adata, lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)), REML = F)

# Plot curve by tertiles of methylation level
cog_adata$gc_predict <-predict(gc_tspat_CpG, na.action = na.exclude)

ggplot(data=cog_adata, aes(x=cog_age, y=tpcom, group=tert_cg180642562, colour=tert_cg180642562)) +
  stat_smooth(aes(y = gc_predict),size=1.2, method = "lm", formula = y ~ x + I(x^2)) +
  labs(x = "Age", y= "Spatial ability", title ="cg18064256 on spatial ability" , color = "Methylation tertile", size=2) +
  scale_color_manual(labels = c("Low", "Middle", "High"), values = c("blue4", "darkcyan", "darkorchid2")) +
  coord_cartesian(xlim=c(55, 90)) +
  theme(plot.title = element_text(size=20, face="bold"),
        axis.title = element_text(size=15, face="bold"),
        axis.text = element_text(size=15),
        legend.justification = c(0, 0), legend.position = c(0, 0),
        legend.title = element_text(size=20, face="bold"),
        legend.text = element_text(size=20)) 
