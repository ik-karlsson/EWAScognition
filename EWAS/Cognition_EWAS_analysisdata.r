##################################################################################################################
# Name: Cognition_EWAS_analysisdata_20191216
#
# Purpose: Create dataframe for analyses, consisting of selected samples
# Study: Cognition EWAS
#
# Created by Ida Karlsson (IK)
# Institute of Gerontology, Jönköping University &
# Department of Medical Epidemiology and Biostatistics (MEB) Karolinska Institutet
#
# Created:  20180504 by Ida Karlsson (IK) (based on script for BMI)
# Updated:	20180504 by IK	- Use full 450K data
#			20181203 by IK  - Add EPIC data
#           20191105 by IK  - Use updated data, IPT 10 added and memory split into
#                             episodic and working memory
#			20191216 by IK  - Converted memory (Thurstone and digit span) to T-scores.
#
#
# Steps:		1. Select samples and create phenotype file
#			    2. Select sites and create file with M-values
#			    3. Save as R data frame
#
##################################################################################################################

.libPaths("/nfs/home/idamat/Tools_and_packages")
library(plyr)

##################################################################################################################
###### 1. Select samples and create phenotype file
##################################################################################################################

load("/nfs/AGE/twins.methylation.data/Epic_blood_SATSA/20171109_1159_pheno_450k_1094.Rdata")
load("/nfs/AGE/twins.methylation.data/Epic_blood_SATSA/20180530_1646_pheno_epic_375.Rdata")
ls()

dim(pheno_450k)
# n=1094x17
dim(pheno_epic)
# n=375x16

head(pheno_450k)
head(pheno_epic)

# Dropping samples from the HARMONY and Parkinson studies
table(pheno_450k$STUDY)
pheno_450k_s <- pheno_450k[ which(pheno_450k$STUDY!='Harmony' & pheno_450k$STUDY!='Parkinson'), ]
dim(pheno_450k_s)
# n=1032x17

# Dropping variables I do not need 
# First convert the sample variable in epic to capital letters..
methvars <- c("Sample", "TWINNR", "PAIRID", "BESTZYG", "STUDY", "SEX", "AGE", "SMOKE")
pheno_450k2 <- pheno_450k_s[methvars]
pheno_epic2 <- pheno_epic[methvars]

# Add chip indicator
pheno_450k2$CHIP <- "450K"
pheno_epic2$CHIP <- "EPIC"

### Stack 450K and EPIC
head(pheno_450k2)
head(pheno_epic2)

meth1 <- rbind(pheno_epic2, pheno_450k2)

##### Select the first available blood sample 
# Sort on twinnr and age
meth1 <- meth1[order(meth1$TWINNR, meth1$AGE),]

# Selecting the first available blood sample
meth2 <- meth1[!duplicated(meth1$TWINNR),]
dim(meth2)
# n=536

########## Importing cognition data from P 
cogdata <- read.table("/cifs/P/Dementia_IK/Cognition/EB_estimates/Data/eb_cognition_all_20191014.txt", header=TRUE, sep="\t")
head(cogdata)

# Add to the meth data
meth3 <- merge(meth2, cogdata, by.x="TWINNR", by.y="twinnr")
dim(meth3)
# n=535x60

# Add indicator for complete pairs
meth4 <- ddply(meth3,.(PAIRID),transform, TWINPAIR = NROW(piece))
table(meth4$TWINPAIR)
# 478 twins belong to complete pair (so 239 complete pairs)

# Saving as pheno_file
pheno_file <- meth4

##################################################################################################################
###### 2. Select samples in the m-value file
##################################################################################################################

# Select m-values present on both 450K and EPIC, with mean beta value difference <0.15
load("/nfs/AGE/twins.methylation.data/Epic_blood_SATSA/arraydif.Rdata")
dim(arraydif)
# n=452453

# Select those with mean diff <0.15
CpG_diff <- arraydif[ which(arraydif$meandif > -0.15 & arraydif$meandif < 0.15 & arraydif$var450),]
dim(CpG_diff)
# n=438643
CpG_list_diff <- as.matrix(rownames(CpG_diff))

# Load beta-values
load("/nfs/AGE/twins.methylation.data/Epic_blood_SATSA/20180128_0046_betas_ct_batch.Rdata")
dim(betas_ct_batch)
# 255356 x 1469

# Select m-values according to the selection list
beta_sites <- betas_ct_batch[rownames(betas_ct_batch) %in% CpG_list_diff, ]
dim(beta_sites)
# n=250816 x 1469

# Select samples according to the selection list
b_values <- beta_sites[,as.character(pheno_file$Sample)]
dim(b_values)
# n=250816 x 535

# Convert to M-values
m_values <- log2(b_values/(1-b_values))

# Read in the Illumina probe annotation (using that from the 450K)
gene_ref_pre <- read.csv(file = "/nfs/AGE/twins.methylation.data/Data/HumanMethylation450_15017482_v1-2.csv",header = TRUE)
rownames(gene_ref_pre) <- gene_ref_pre$IlmnID
head(gene_ref_pre)
# Select sites included in the mval_sites data
gene_ref <- merge(m_values, gene_ref_pre, by="row.names", sort=F)
# Assign the CpG ID as rowname, and keep only gene name
rownames(gene_ref) <- gene_ref$Row.names
gene_ref_var <- c("UCSC_RefGene_Name", "CHR", "MAPINFO")
gene_ref <- gene_ref[gene_ref_var]
head(gene_ref)
dim(gene_ref)
# 250816 x 3

##################################################################################################################
##### 3. Save as R data frame

# Some checks:
m_check <- m_values[,1:5]
p_check <- pheno_file[c(1:5)]

head(p_check)
head(m_check)
head(gene_ref)

tail(m_check)
tail(gene_ref)

tail(m_values)
tail(p_check)


save(m_values, pheno_file, gene_ref, file="/nfs/home/idamat/EWAS_cognition/Data/Cognition_EWAS_analysisdata_20191216.Rdata")

##################################################################################################################
##### END OF FILE
##################################################################################################################
