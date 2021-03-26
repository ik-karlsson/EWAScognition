##################################################################################################################
# Name: Cognition_EWAS_twincorr

# Purpose: Look at twin-pair correlation of CpG sites
# Study: EWAS of cognition

# Created by Ida Karlsson (IK)
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20210310 by Ida Karlsson (IK) 
# Updated:	

# Input: P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20191216.txt
#        P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tmemod_20191216.txt

# Output: P:/Dementia_IK/Cognition/EWAS/Output/Reviewer_comments/Cognition_EWAS_twincorr_20210310.txt

##################################################################################################################

##################################################################################################################
########## Load libraries
.libPaths("Z:/Programs/R/Packages")
library(stringr)

adata1 <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tspat_20191216.txt", header=T)
adata1 <- adata1[c("TWINNR", "PAIRID", "BESTZYG", "cg04549090", "cg18064256")]

adata2 <- read.table("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_mval_65_tmemod_20191216.txt", header=T)
adata2 <- adata2[c("TWINNR", "cg09988380", "cg25651129", "cg08011941")]

adata3 <- merge(adata1, adata2, by="TWINNR")

# Create TVAB
adata3$TVAB <- substr(adata3$TWINNR, nchar(adata3$TWINNR), nchar(adata3$TWINNR))

# Subset into tvab 1 and 2
# Drop variables first
tvab1 <- adata3[which(adata3$TVAB==1),]
tvab2 <- adata3[which(adata3$TVAB==2),]

# Merge by pairid
pairs <- merge(tvab1, tvab2, by="PAIRID")

# Subset by zygosity
MZpairs <- pairs[which(pairs$BESTZYG.x==1),]
DZpairs <- pairs[which(pairs$BESTZYG.x==2),]

##### Create table to store results
corr_results <- matrix(NA, nrow = 6, ncol = 4)

# Assign some general columns
corr_results[1, 1:4] = c("CpG", "rMZ (95% CI)", "rDZ (95% CI)", "Falconer heritability")

corr_results[c(2:6), 1] = colnames(adata3)[4:8]

##### Loop over Cpg sites and output correlations to results file
for (j in 1:5){
  
  rMZ <- cor.test(MZpairs[,(j+3)], MZpairs[,(j+11)], method = "pearson")
  rDZ <- cor.test(DZpairs[,(j+3)], DZpairs[,(j+11)], method = "pearson")
  
  est_MZ = formatC(rMZ$estimate, digits = 2, format = "f")
  est_DZ = formatC(rDZ$estimate, digits = 2, format = "f")
  
  ci_MZ = formatC(rMZ$conf.int, digits = 2, format = "f")
  ci_DZ = formatC(rDZ$conf.int, digits = 2, format = "f")
  
  corr_results[j+1, 2] = str_trim(paste(est_MZ, " (", ci_MZ[1], "-", ci_MZ[2], ")", sep=""))
  corr_results[j+1, 3] = str_trim(paste(est_DZ, " (", ci_DZ[1], "-", ci_DZ[2], ")", sep=""))
  
  corr_results[j+1, 4] = formatC(2*(rMZ$estimate-rDZ$estimate), digits = 2, format = "f")
  
}

##### Output results

write.table(corr_results, "P:/Dementia_IK/Cognition/EWAS/Output/Reviewer_comments/Cognition_EWAS_twincorr_20210310.txt", col.names = F, row.names = F, quote = F, sep="\t")

##################################################################################################################
##################################################################################################################
##################################################################################################################


