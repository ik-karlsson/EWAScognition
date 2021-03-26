##################################################################################################################
# Name: Extract_CpG_Mvalues_BMI

# Purpose: Extract M-values for significant CpGs in the cognition EWAS
#
# Study: EWAS of cognition
#
# Created by Ida Karlsson (IK)
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet, Stockholm, Sweden
#
# Created:   20210324 by Ida Karlsson (IK)
# Updated:   
#
# Based on code from BMI methylation study. Original code from Yunzhang Wang
#
##################################################################################################################

##### OBTAIN M-VALUES
load("/nfs/home/idamat/EWAS_cognition/Data/Cognition_EWAS_analysisdata_20191216.Rdata")
dim(m_values)
# 250816  x  535

# Select the 5 CpG sites
CpG_selection <- c("cg04549090", "cg18064256", "cg09988380", "cg25651129", "cg08011941" )
m_selection <- as.data.frame(m_values[rownames(m_values) %in% CpG_selection, ])
dim(m_selection)
# 5 x 535

# Transpose to long format
t_m_selection <- t(m_selection)
head(t_m_selection)
dim(t_m_selection)
# 535 x 5

# Add twinnr from the pheno file
twinnr <- pheno_file[,1:2]
m_twinnr <- merge(twinnr, t_m_selection, by.x="Sample", by.y="row.names")
head(m_twinnr)
m_twinnr <- m_twinnr[-1]

# Save as text file
write.table(m_twinnr, "/nfs/home/idamat/EWAS_cognition/Data/Extract_CpG_Mvalues_20210324.txt", row.names = F, quote=F)

### END OF FILE
