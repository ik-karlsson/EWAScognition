##################################################################################################################
# Name: Cognition_EWAS_CpG_followup_20181217

# Purpose: Select all nominally significant CpGs from the EWAS runs for follow-up analyses
# Study: Cognition EWAS

# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet

# Created:  20181217 by Ida Karlsson (IK) (based on script for BMI)
# Updated:	20191014 by IK - New data incl IPT10
#							 Split memory into Thurstone and digit span.  
#           20191216 by IK - T-scores for memory, standardized M-values, and weighted regression
#							 Add pair means and deviations from pair means for between-within analyses
#							 Add chromosome and position

# Input: Output files from Cognition_EWAS_analyses_20191216

##################################################################################################################

### Load the m-value data
# Read the R dataframe contianing m-values(m_values)
load("/nfs/home/idamat/EWAS_cognition/Data/Cognition_EWAS_analysisdata_20191216.Rdata")

# Select relevant columns from the pheno file (will be merged on each set of CpGs)
pheno_clean <- pheno_file[c(1:9)]
dim(pheno_clean)
# n=535x9
	
##################################################################################################################
##### Loop over verbal ability, spatial ability, processing speed, memory/thurstone, and general cognition
##### NOTE! memory/digit span below as it was a linear model
##### For each, the EB estimate for intercept, linear slope, and quadratic slope

# Generate list of outcomes
varlist <- c("70_tverb", "65_tsped", "65_tmemot", "65_tspat", "65_tpcom")

# Loop over list
for(i in 1:5) {
  outcome <- varlist[i]	
  
  ### Intercept
  # Specify the input file
  infile_icpt <- paste("/nfs/home/idamat/EWAS_cognition/Results/qI",outcome,"_adj_top_20191216.txt",sep="")
  results_icpt <- read.table(infile_icpt, header=T, sep=" ")
  results_icpt$EBpiece <- "Intercept"
  # Select CpGs with p<10-5
  sign_icpt <- results_icpt[ which(results_icpt$P_value < 10^-5), ]
  
  ### Linear slope
  # Specify the input file
  infile_linear <- paste("/nfs/home/idamat/EWAS_cognition/Results/qL",outcome,"_adj_top_20191216.txt",sep="")
  results_linear <- read.table(infile_linear, header=T, sep=" ")
  results_linear$EBpiece <- "Linear slope"
  # Select CpGs with p<10-5
  sign_linear <- results_linear[ which(results_linear$P_value < 10^-5), ]
  
  ### Quadratic slope
  # Specify the input file
  infile_quad <- paste("/nfs/home/idamat/EWAS_cognition/Results/qQ",outcome,"_adj_top_20191216.txt",sep="")
  results_quad <- read.table(infile_quad, header=T, sep=" ")
  results_quad$EBpiece <- "Quadratic slope"
  # Select CpGs with p<10-5
  sign_quad <- results_quad[ which(results_quad$P_value < 10^-5), ]

  # Stack the files
  sign_all <- rbind(sign_icpt, sign_linear, sign_quad)
  
  # Add chromosome and position
  sign_all_cp <- merge(gene_ref, sign_all, by.x="row.names", by.y="CpG")
  colnames(sign_all_cp)[colnames(sign_all_cp)=="Row.names"] <- "CpG"
  sign_all_cp <- sign_all_cp[c(-2)]
  
  # Output
  output <- paste("/nfs/home/idamat/EWAS_cognition/Results/CpG_followup_table_",outcome,"_20191216.txt",sep="")
  write.table(sign_all_cp, output, row.names = F)
  
  ### Select the relevant CpGs from the mval file
  CpG_list <- sign_all$CpG
  sign_mval <- m_values[rownames(m_values) %in% CpG_list, ]
  
  # Transpose to long format and standardize
  t_sign_mval <- t(sign_mval)
  for (j in 1:ncol(t_sign_mval)) {
	t_sign_mval[,j] <- scale(t_sign_mval[,j])
  }
  
  # Add EB estimates for the relevant outcome
  icpt <- paste("qI",outcome,sep="")
  Lslope <- paste("qL",outcome,sep="")
  Qslope <- paste("qQ",outcome,sep="")
  icpt_w <- paste("qI",outcome,"se_inv",sep="")
  Lslope_w <- paste("qL",outcome,"se_inv",sep="")
  Qslope_w <- paste("qQ",outcome,"se_inv",sep="")

  add_pheno <- pheno_file[c("TWINNR", icpt, icpt_w, Lslope, Lslope_w, Qslope, Qslope_w)]
  
  # Merge on to the pheno clean file
  pheno_out <- merge(pheno_clean, add_pheno, , by="TWINNR")
  
  # And the M-values
  mval_pheno <- merge(pheno_out, t_sign_mval, , by.x="Sample", by.y="row.names")
  mval_pheno <- mval_pheno[order(mval_pheno$TWINNR),]
  
  # Generate file with pair means and deviations from pair means for between-within analyses 
  nCpG=nrow(sign_mval)
  # Generate pair mean
  for (k in 1:nCpG){
    mval_pheno[,k+15+nCpG] <- ave(mval_pheno[,k+15], as.numeric(mval_pheno$PAIRID), FUN=mean)
    colnames(mval_pheno)[k+15+nCpG] <- paste("pairm_", colnames(mval_pheno)[k+15], sep = "")
  }
  # Create departure from mean pair m-value for each twin
  for (l in 1:nCpG){
    mval_pheno[,l+15+2*nCpG] <- mval_pheno[,l+15]-mval_pheno[,l+15+nCpG]
    colnames(mval_pheno)[l+15+2*nCpG] <- paste("depm_", colnames(mval_pheno)[l+15], sep = "")
  }
  # Output
  output_mval <- paste("/nfs/home/idamat/EWAS_cognition/Results/CpG_followup_mval_",outcome,"_20191216.txt",sep="")
  write.table(mval_pheno, output_mval, na=".", sep="\t", quote=F, row.names=F)
}

##################################################################################################################
##### Memory/digit span 
##### EB estimate for intercept and linear slope 

### Intercept
# Specify the input file
results_icpt <- read.table("/nfs/home/idamat/EWAS_cognition/Results/qI65_tmemod_adj_top_20191216.txt", header=T, sep=" ")
results_icpt$EBpiece <- "Intercept"
# Select CpGs with p<10-5
sign_icpt <- results_icpt[ which(results_icpt$P_value < 10^-5), ]

### Linear slope
# Specify the input file
results_linear <- read.table("/nfs/home/idamat/EWAS_cognition/Results/qL65_tmemod_adj_top_20191216.txt", header=T, sep=" ")
results_linear$EBpiece <- "Linear slope"
# Select CpGs with p<10-5
sign_linear <- results_linear[ which(results_linear$P_value < 10^-5), ]

# Stack the files
sign_all <- rbind(sign_icpt, sign_linear)

# Add chromosome and position
  sign_all_cp <- merge(gene_ref, sign_all, by.x="row.names", by.y="CpG")
  colnames(sign_all_cp)[colnames(sign_all_cp)=="Row.names"] <- "CpG"
  sign_all_cp <- sign_all_cp[c(-2)]
  
# Output
output <- paste("/nfs/home/idamat/EWAS_cognition/Results/CpG_followup_table_65_tmemod_20191216.txt",sep="")
write.table(sign_all_cp, output, row.names = F)

### Select the relevant CpGs from the mval file
CpG_list <- sign_all$CpG
sign_mval <- m_values[rownames(m_values) %in% CpG_list, ]

# Transpose to long format and standardize
t_sign_mval <- t(sign_mval)
for (j in 1:ncol(t_sign_mval)) {
	t_sign_mval[,j] <- scale(t_sign_mval[,j])
}
  
add_pheno <- pheno_file[c("TWINNR", "qI65_tmemod", "qI65_tmemodse_inv", "qL65_tmemod", "qL65_tmemodse_inv")]
  
# Merge on to the pheno clean file
pheno_out <- merge(pheno_clean, add_pheno, , by="TWINNR")
  
# And the M-values
mval_pheno <- merge(pheno_out, t_sign_mval, , by.x="Sample", by.y="row.names")
mval_pheno <- mval_pheno[order(mval_pheno$TWINNR),]
  
# Generate file with pair means and deviations from pair means for between-within analyses 
nCpG=nrow(sign_mval)
# Generate pair mean
for (k in 1:nCpG){
  mval_pheno[,k+13+nCpG] <- ave(mval_pheno[,k+13], as.numeric(mval_pheno$PAIRID), FUN=mean)
  colnames(mval_pheno)[k+13+nCpG] <- paste("pairm_", colnames(mval_pheno)[k+13], sep = "")
}
# Create departure from mean pair m-value for each twin
for (l in 1:nCpG){
  mval_pheno[,l+13+2*nCpG] <- mval_pheno[,l+13]-mval_pheno[,l+13+nCpG]
  colnames(mval_pheno)[l+13+2*nCpG] <- paste("depm_", colnames(mval_pheno)[l+13], sep = "")
}

# Output
write.table(mval_pheno, "/nfs/home/idamat/EWAS_cognition/Results/CpG_followup_mval_65_tmemod_20191216.txt", na=".", sep="\t", quote=F, row.names=F)

##################################################################################################################
