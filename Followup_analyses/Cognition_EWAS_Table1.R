##################################################################################################################
# Name: Cognition_EWAS_Table1_20200113

# Purpose: Generate Table 1 for MS
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20200113 by Ida Karlsson (IK) (based on script for BMI)
# Updated:	20200820 by IK - Scale linear and quadratic estimates to 10-year change

##################################################################################################################
########## Load libraries
.libPaths("Z:/Programs/R/Packages")
library(dplyr)

# Generate files for significant and suggestive sites
# Specify header
addrow1 <- c("Outcome"," ","CpG information"," "," "," ","Main EWAS analysis"," "," ","Between-pair effect"," "," ","Within-pair effect"," "," ")
addrow2 <- c("Cognitive domain",  "Trajectory feature", "CpG", "Gene", "Chromosome", "Position","Beta", "SE", "P-value", "Beta", "SE", "P-value", "Beta", "SE", "P-value")
results_sign <- as.matrix(rbind(addrow1, addrow2))
results_sugg <- rbind(addrow1, addrow2)

# Generate list of outcomes
varlist <- c("65_tsped", "70_tverb", "65_tspat", "65_tmemot", "65_tmemod", "65_tpcom")
out_list <- c("Processing speed", "Verbal ability", "Spatial ability", "Episodic memory", "Working memory", "General cognitive ability")

# Loop over list
for(i in 1:6) {
  outcome <- varlist[i]	
  domain <- substring(outcome, 4, nchar(outcome))
  # Specify the input files
  # Main results
  infile_main <- paste("P:/Dementia_IK/Cognition/EWAS/Output/CpG_followup_table_",outcome,"_20191216.txt",sep="")
  results_main <- read.table(infile_main, header=T, sep=" ", stringsAsFactors=F)
  # Split the gene ref into just one (no repetitions)
  results_main$Gene <- sub(';.*$','', results_main$Gene_ref)
  
  # Between-within
  infile_bw <- paste("P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_BetweenWithin_20200106_",domain,".txt",sep="")
  results_bw <- read.table(infile_bw, header=T, sep="\t")
  
  # Merge
  results_all <- merge(results_main, results_bw, by.x=c("CpG", "EBpiece"), by.y=c("CpG", "outcome"))

  # Scale the linear and quadratic estimates to represent 10-year change
  for(j in c(6,7,11,12,14,15)){
    results_all[,j] <- ifelse((results_all$EBpiece =="Intercept"), results_all[,j],
                              ifelse((results_all$EBpiece =="Linear slope"), results_all[,j]*10,
                                     ifelse((results_all$EBpiece =="Quadratic slope"), results_all[,j]*100, NA)))
  }
  
    
  # Sort by EB-piece and chromosome, and position
  results_all <- results_all[order(results_all$EBpiece, results_all$CHR, results_all$MAPINFO),]
  
  # Add cognitive domain
  results_all$outcome <- out_list[i]
  
  # Order columns
  results_all_c <- results_all[c(17,2,1,10,3,4,6:8,11:16)]

  # Split into significant and suggestive
  add_sign <- results_all_c[which(results_all_c$P_value<2.4e-7),]
  add_sugg <- results_all_c[which(results_all_c$P_value>=2.4e-7),]
  
  # Significant sites (if any)
  if (nrow(add_sign)>=1){
    # Deal with formatting
    # Betas and SEs
    for (j in c(7,8,10,11,13,14)) {
      add_sign[,j] <- formatC(add_sign[,j], digits = 2, format = "f")
    }
    # P-values
    for (k in c(9,12,15)) {
      add_sign[,k] <- format.pval(add_sign[,k], sci=-2, digits=3)
    }
    # Convert to matrix and rbind with the output file
    add_sign_m <- as.matrix(add_sign)
    results_sign <- rbind(results_sign, add_sign_m)
  }

  # Same for suggestive sites
  # Deal with formatting
  # Betas and SEs
  for (l in c(7,8,10,11,13,14)) {
    add_sugg[,l] <- formatC(add_sugg[,l], digits = 2, format = "f")
  }
  # P-values
  for (m in c(9,12,15)) {
    add_sugg[,m] <- format.pval(add_sugg[,m], sci=-2, digits=3)
  }
  
  # Convert to matrix and rbind with the output file
  add_sugg_m <- as.matrix(add_sugg)
  results_sugg <- rbind(results_sugg, add_sugg_m)
  
}  

# Output as text file
write.table(results_sign, "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_Table1_sign_20200820.txt", row.names=F, col.names=F, sep="\t", quote=F)
write.table(results_sugg, "P:/Dementia_IK/Cognition/EWAS/Results/Cognition_EWAS_Table1_sugg_20200820.txt", row.names=F, col.names=F, sep="\t", quote=F)



