## ----------------------------------------------
## Filename: preprocess_norm.R
## Study: EPIC_blood
## Author: Yunzhang Wang
## Date: 2017-09-28
## Updated:
## Purpose: The function to normalize EPIC chip data
## Notes: minfi methods
## -----------------------------------------------
## Data used:
## * QCed RGset
## -----------------------------------------------
## Output:
## * Function normalized Mset
## -----------------------------------------------
## OP: R 3.4.1
## -----------------------------------------------*/

preprocess_norm <- function(
  RGset_qc,
  bgcorrect=T,
  savedata=F,
  filename="Mset_funnorm.Rdata",
  norm.method="funnorm",
  .out.dir="/proj/b2016305/INBOX/PE-1124/PE-1124_170823_ResultReport/output"

){

require(minfi)
source("~/epic_aging/src/preload/global_functions.R")

########## Functional normalization

if(norm.method=="funnorm"){
  if(bgcorrect)  TerminalLogger("Start background correction.")
  TerminalLogger("Start minfi functional normalization.")
  Mset_norm <- preprocessFunnorm(RGset_qc, nPCs=2, , bgCorr = bgcorrect, dyeCorr = TRUE, verbose = TRUE) # function normalization
}

if(norm.method=="dasen") {
  require(wateRmelon)
  if(bgcorrect) {
    TerminalLogger("Start background correction.")
    RGset_qc <- preprocessNoob(RGset_qc)
  }
  TerminalLogger("Start wateRmelon dasen normalization.")
  Mset_norm <- dasen(RGset_qc)

}
if(savedata) save(Mset_norm, file=file.path(.out.dir, filename))

TerminalLogger("Finish normalization.")


return(Mset_norm)
}
