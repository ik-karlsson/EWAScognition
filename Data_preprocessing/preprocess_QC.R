## ----------------------------------------------
## Filename: preprocess_QC.R
## Study: EPIC_blood
## Author: Yunzhang Wang
## Date: 2017-09-28
## Updated:
## Purpose: The function to remove bad samples and probes
## Notes: minfi methods
## -----------------------------------------------
## Data used:
## * Raw RGset
## -----------------------------------------------
## Output:
## * QCed RGset
## -----------------------------------------------
## OP: R 3.4.1
## -----------------------------------------------*/


preprocess_QC <- function(
  ifloadrgset = T,
  RGset_raw = NULL,
  dpval = 0.01,
  rmsexProbe=T, # remove probes on sex chromosomes or not
  rmsnpProbe=T, # remove probes with a snp or not
  rmoutliers=F, # remove outliers in PC1
  keepcpg = NULL,
  rmsample = NULL, # 201496850093_R07C01 201496860021_R03C01 201503670171_R01C01 low quality
                   # 201496850061_R06C01 201496860035_R01C01 bad distribution
  keepsample = NULL,
  only450k = F,
  saverawbeta= F,
  .out.dir = "/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/output",
  savedata=F,
  filename="RGset_qc.Rdata") {

require(minfi)
require(IlluminaHumanMethylationEPICmanifest)
require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
source("/home/yunzhang/epic_aging/src/preload/global_functions.R")

raw.dir <- "/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/raw"
out.dir <- "/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/output"

########## load data
if(ifloadrgset){
  load(file.path(out.dir, "RGset_raw.Rdata"))
  RGset_raw <- RGset
}

betas_raw <- getBeta(RGset_raw)

if(saverawbeta){
  save(betas_raw, file = file.path(.out.dir, TimeStampMe("betas_raw.Rdata")))
}

pval <- detectionP(RGset_raw)   ### get detectionP-values from raw rgset
gc()

########## QC
### remove samples
samples <- rep(F, ncol(pval))
names(samples) <- colnames(pval)
keepsample <- keepsample[keepsample %in% names(samples)]

if(!only450k){
  ss <- read.delim(file.path(out.dir, "sample_sheet.txt"))
  rownames(ss) <- ss$Sample


  bloodsamples <- as.vector(ss[grep("SATSA", ss$Sample.ID), "Sample"])
  samples[bloodsamples[bloodsamples %in% names(samples)] ] <- T

  if(rmoutliers) {
    samples["201496860021_R07C01"] <- F
    samples["201533510019_R04C01"] <- F
    samples["201533490006_R08C01"] <- F
  }

  samples[rmsample] <- F
}
samples[keepsample] <- T

TerminalLogger(paste0("Remove ",length(samples[!samples])," samples, ", length(samples[samples]), " samples left." ))

### remove probes
badProbe1 <- apply(pval[,samples] > dpval, 1, any)    # Determine probes having p-val over 0.05 or 0.01 in any sample
                                                      # 165196 CpGs (high), 297359 CpGs (intermediate)
badProbe1 <- names(badProbe1[badProbe1])

badProbe2 <- apply(is.na(betas_raw[,samples]), 1, any)    # Determine probes having missing values in any sample
badProbe2 <- names(badProbe2[badProbe2])

badProbe <- unique(c(badProbe1, badProbe2))
TerminalLogger(paste0("Remove ",length(badProbe), " probes, including ", length(badProbe1)," low-quality probes and ", length(badProbe2), " missing-value probes." ))

### Remove probe-overlapping SNP
if(rmsnpProbe){
  snp <- getSnpInfo(RGset_raw)  # get info of probe-overlapping SNP
  snpProbe <- rownames(subset(as.data.frame(snp), !is.na(CpG_maf) | !is.na(SBE_maf) | !is.na(Probe_maf))) # 181462 CpGs
  TerminalLogger(paste0("Remove ",length(snpProbe)," SNP-overlapping probes."))
} else {
  snpProbe <- NULL
}
# snpProbe.1 <- rownames(subset(as.data.frame(snp), !is.na(CpG_maf) | !is.na(SBE_maf) | (!is.na(Probe_maf) & Probe_maf > 0.05) ))

### Remove sex probes
if(rmsexProbe){
  sexProbe <- rownames(subset(as.data.frame(getLocations(RGset_raw)), seqnames %in% c("chrX", "chrY"))) # 19681 CpGs
} else{
  sexProbe <- NULL
}

### perform QC
rmProbe <- unique(c(badProbe, snpProbe, sexProbe))
rmProbe <- rmProbe[!rmProbe %in% keepcpg]
gc()

RGset_qc <- subsetByLoci(RGset_raw[, samples], excludeLoci = rmProbe,keepControls = TRUE, keepSnps = TRUE)  # remove probes&samples from RGset

RGset_qc

if(savedata) save(RGset_qc, file=file.path(.out.dir, filename))

return(RGset_qc)
gc()

}
# if(rmsnpProbe) write.table(as.data.frame(snp), file=file.path(out.dir, "SNP-overlapping.txt"), sep="\t", quote=F, row.names=F)
