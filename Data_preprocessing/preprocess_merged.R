## ----------------------------------------------
## Filename: preprocess_merged.R
## Study: EPIC_blood
## Author: Yunzhang Wang
## Date: 2017-09-28
## Updated:
## Purpose: To preprocess merged EPIC and 450k chip data
## Notes: depends on preprocess_QC.R, preprocess_norm.R preprocess_cellbatch.R
## -----------------------------------------------
## Data used:
## * raw RGset
## -----------------------------------------------
## Output:
## * beta-values
## -----------------------------------------------
## OP: R 3.4.1
## -----------------------------------------------*/

########## settings
wd.dir <- "/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/output/combine/sample1469_p005_dasen"
if(!dir.exists(wd.dir)) dir.create(wd.dir)
out.dir <- "/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/output"
in.dir <- "/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/output/combine"

require(lumi)
require(minfi)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

source("/home/yunzhang/epic_blood/src/preprocess_QC.R")
source("/home/yunzhang/epic_blood/src/preprocess_norm.R")
source("/home/yunzhang/epic_blood/src/preprocess_cellbatch.R")

phe <- load_data(uppmax=T, pheno=T, greedy.cut=F)
phe <- phe$pheno
phe$name <- paste0(phe$SLIDE, "_", phe$POSITION)

########## load data
load(file.path(out.dir, "RGset_raw.Rdata"))
RGset_epic <- RGset

load(file.path(in.dir, "RGset_raw.Rdata"))
RGset_450 <- RGset

RGset_cb <- combineArrays(RGset_epic, RGset_450)

rm(RGset_epic)
rm(RGset_450)
gc()

########## QC
qced <- preprocess_QC(
  ifloadrgset = F,
  RGset_raw = RGset_cb,
  dpval = 0.05,
  rmsexProbe=T, # remove probes on sex chromosomes or not
  rmsnpProbe=T, # remove probes with a snp or not
  rmoutliers=T, # remove outliers in PC1
  keepcpg = NULL,
  rmsample = c("201496850061_R06C01", "201496860035_R01C01"), # 201496850093_R07C01 201496860021_R03C01 201503670171_R01C01 low quality
                   # 201496850061_R06C01 201496860035_R01C01 bad distribution
  keepsample = as.vector(phe$name),
  # keepsample = colnames(RGset_cb)[6:10],
  saverawbeta= T,
  savedata=T,
  .out.dir=wd.dir)

cellcounts <- estimate_cellcounts(qced, .out.dir=wd.dir)
dim(cellcounts)
save(cellcounts, file=file.path(wd.dir, TimeStampMe("cellcounts.Rdata")))

normed <- preprocess_norm(
  RGset_qc=qced,
  bgcorrect=T,
  savedata=F,
  norm.method="dasen",
  .out.dir=wd.dir
)

bt <- normed@colData@listData$Slide

betas <- getBeta(normed)  #  CpG,  samples
save(betas, file=file.path(wd.dir, TimeStampMe("betas_norm.Rdata")))

rm(qced)
rm(normed)
gc()

summary(colnames(betas) == rownames(cellcounts))

betas_ct <- adjust_cellcounts(.Beta=t(betas), celltypes=cellcounts)
betas_ct <- t(betas_ct)
mval_ct <- beta2m(betas_ct)

save(betas_ct, file=file.path(wd.dir, TimeStampMe("betas_ct.Rdata")))

rm(betas_ct)
gc()

mval_ct_batch <- adjust_batch(edata=mval_ct, batches=bt, .out.dir=wd.dir, parp=T, name="_ct_batch")
betas_ct_batch <- m2beta(mval_ct_batch)

save(mval_ct_batch, file=file.path(wd.dir, TimeStampMe("mval_ct_batch.Rdata")))
save(betas_ct_batch, file=file.path(wd.dir, TimeStampMe("betas_ct_batch.Rdata")))

rm(mval_ct_batch, betas_ct_batch)
gc()

mval <- beta2m(betas)
mval_batch <- adjust_batch(edata=mval, batches=bt, .out.dir=wd.dir, parp=T, name="_batch")
betas_batch <- m2beta(mval_batch)
save(mval_batch, file=file.path(wd.dir, TimeStampMe("mval_batch.Rdata")))
save(betas_batch, file=file.path(wd.dir, TimeStampMe("betas_batch.Rdata")))
