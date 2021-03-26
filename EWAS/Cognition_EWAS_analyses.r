##################################################################################################################
# Name: Cognition_EWAS_analyses_20191216
#
# Purpose: Run EWAS analyses on Cognitive domains
# Study: Cognition EWAS
#
# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet
#
# Created:  20180508 by Ida Karlsson (IK) (based on script for BMI)
# Updated:	20180504 by IK	- Use full 450K data, run adjusted models only
#			20181203 by IK  - Use 450K + EPIC data. Add models with linear slope
#			20181207 by IK  - Create one big loop over all outcomes
#           20191105 by IK  - Use updated data, IPT 10 added and memory split into
#                             episodic and working memory
#			20191216 by IK  - Converted memory (Thurstone and digit span) to T-scores.
#							  Include weights to account for number of measurements
#
# Using EB estimates for cognitive domains, modelling the intercept, linear slope, and quadratic slope
#
# Based on code for EWAS of dementia
#
# Input: /nfs/home/idamat/EWAS_cognition/Data/Cognition_EWAS_analysisdata_20191216.Rdata
#
##################################################################################################################

# Load required libraries
.libPaths("/nfs/home/idamat/Tools_and_packages")
library(Formula)
library(multiwayvcov)
library(lmtest)

# Read the R dataframe contianing phenotype information (pheno_file) and m-values(m_values)
# NOTE! Phenotype data is in long format (ID in rows), while the m-value data is in wide format (ID in columns)
load("/nfs/home/idamat/EWAS_cognition/Data/Cognition_EWAS_analysisdata_20191216.Rdata")

##################################################################################################################
##### Loop over all outcomes
##### Verbal ability, spatial ability, processing speed, episodic memory, working memory, and general cognition
##### For each, the EB estimate for intercept, linear slope, and quadratic slope

# Generate list of outcomes
varlist <- c("qI70_tverb", "qL70_tverb", "qQ70_tverb", "qI65_tsped", "qL65_tsped", "qQ65_tsped", 
			"qI65_tmemot", "qL65_tmemot", "qQ65_tmemot", "qI65_tmemod", "qL65_tmemod", 
			"qI65_tspat", "qL65_tspat", "qQ65_tspat", "qI65_tpcom", "qL65_tpcom", "qQ65_tpcom")

# Loop over list
for(i in 1:18) {
  outcome <- varlist[i]	
  weight <- paste(outcome, "se_inv", sep="")

# Create a matrix of the values and number of columns to stored
lr_results_adj <- matrix(NA, nrow = nrow(m_values), ncol = 4)
colnames(lr_results_adj) <- c("Estimate", "SE", "P_value", "P_value_adj")
lr_results_adj = data.frame(CpG      = rownames(m_values), 
                              Gene_ref = gene_ref$UCSC_RefGene_Name, lr_results_adj)

## Linear regression model of M-values on the cognitive outcome variables
### Adjusted model with age, sex, and smoking as covariates (robust standard errors to account for relatedness among twins)

for (j in 1:nrow(m_values))
{
	# Standardized M-values
	X <- (m_values[j, ]	- mean(m_values[j,]))/sd(m_values[j,]) 
	lr_CpG_adj <- lm(get(outcome) ~ X + AGE + SEX + SMOKE + CHIP, data=pheno_file, weights=get(weight))
	lr_CpG_adj.vcovCL<-cluster.vcov(lr_CpG_adj, as.vector(pheno_file$PAIRID)) 
	coefficients <- coeftest(lr_CpG_adj, lr_CpG_adj.vcovCL)
	# Store the results in the matrix defined above
	lr_results_adj[j, 3:5] = coefficients[2,c(1,2,4)]
}				

# Adjust p-values using FDR and add the results to the results matrix
lr_results_adj[,6] = p.adjust(lr_results_adj[,5], method="BH")		

# Sort by p-value
lr_results_adj <- lr_results_adj[order(lr_results_adj[, "P_value"]) ,]

# Select first 100 rows
lr_results_adj_100 <- lr_results_adj[1:100,]

# Save as text file
# Full data
output_full <- paste("/nfs/home/idamat/EWAS_cognition/Results/",outcome,"_adj_full_20191216.txt",sep="")
write.table(lr_results_adj, output_full, row.names = F)
# Topp 100
output_top <- paste("/nfs/home/idamat/EWAS_cognition/Results/",outcome,"_adj_top_20191216.txt",sep="")
write.table(lr_results_adj_100, output_top, row.names = F)
}

q()
n

##################################################################################################################
##### END OF FILE
##################################################################################################################
