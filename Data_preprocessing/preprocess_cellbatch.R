source("~/epic_aging/src/preload/global_functions.R")

########## estimate cellular compositions
estimate_cellcounts <- function(.RGset, .out.dir){
  require(minfi)
  try(.RGset@colData@listData$filenames <- NULL) # Remove extra column to avoid error

  pdf(file.path(.out.dir, "celltypes.pdf"))
  cellcounts <- estimateCellCounts(.RGset, cellTypes=c("Bcell","CD4T","CD8T","Eos","Gran","Mono","NK"), meanPlot=T) # May try different cell types to avoid error
  dev.off()

  return(cellcounts)

}

########## adjust cellular compositions
adjust_cellcounts <- function(.Beta, celltypes){
  TerminalLogger("Correct for celltype composition!")
  ## Get residuals
  celltypes.matrix <- as.matrix(celltypes)  # matrix faster, remove name-row
  cell.residuals <- lm(.Beta ~ celltypes.matrix)$residuals  # sample x variable
###  Adjust Beta
  TerminalLogger("Create adjusted betas")
  ## Adjust Beta! Add residuals of each regression model to the mean methylation
  ## value of each probe (mean across all samples) to obtain the adjusted
  ## methylation data.
  mean.matrix <- matrix(colMeans(.Beta),  # probe mean vector
                        nrow = nrow(cell.residuals), ncol = ncol(cell.residuals),
                        byrow = TRUE)     # nrow x probes
  Beta_adj <- cell.residuals + mean.matrix  # sample x probe
  ## Fix betas outside of allowed range: ]0,1[
  Beta_adj <- apply(Beta_adj, 2, function(probe.col){
      old.col <- probe.col
      ## for beta belonging to a probe
      ## if beta > 1: Set to largest value < 1 in this probe set of betas
      probe.col[probe.col > 1] <- max(subset(old.col, old.col < 1))
      ## if beta < 0: Set to smallest value > 0 in this probe set of betas
      probe.col[probe.col < 0] <- min(subset(old.col, old.col > 0))
      return(probe.col)
  })
  ## Convert to M-values for Combat
  if((all(Beta_adj < 1)) & (all(Beta_adj > 0))) {
    return(Beta_adj)
  } else {
      stop("Invalid beta-values: !(0 < beta < 1)")
  }
}

## adjust cellcounts on mval
adjust_cellcounts_mval <- function(.mval, celltypes){
  TerminalLogger("Correct for celltype composition!")
  ## Get residuals
  celltypes.matrix <- as.matrix(celltypes)  # matrix faster, remove name-row
  cell.residuals <- lm(.mval ~ celltypes.matrix)$residuals  # sample x variable
###  Adjust mval
  TerminalLogger("Create adjusted betas")
  ## Adjust mval! Add residuals of each regression model to the mean methylation
  ## value of each probe (mean across all samples) to obtain the adjusted
  ## methylation data.
  mean.matrix <- matrix(colMeans(.mval),  # probe mean vector
                        nrow = nrow(cell.residuals), ncol = ncol(cell.residuals),
                        byrow = TRUE)     # nrow x probes
  mval_adj <- cell.residuals + mean.matrix  # sample x probe
}

########## adjust batch effects
adjust_batch <- function(edata , batches, .out.dir, parp=T, name=NULL, meanonly=F){
  require(sva)

  modcombat <- model.matrix(~1, data=factor(batches))

  if(any(table(batches) < 2)){
    message("Cannot correct batch effect on SLIDE due to limited samples")
    return(edata)

  } else{
    TerminalLogger("Start correcting batch effect.")
    png(file = file.path(.out.dir, TimeStampMe(paste0("Combat_plot",name,".png"))))
    mval.preprocessed <-  # probe x sample
    sva::ComBat(dat = edata, batch = batches, mod = modcombat,
      ## par.prior = F performs nonparametric empirical
      ## Bayesian adjustments (slower). prior.plots shows
      ## estimated difference between par.prior == T OR F.
      par.prior = parp,
      mean.only=meanonly,
      prior.plots=T)
    dev.off()
    TerminalLogger("Finish correcting batch effect.")
    return(mval.preprocessed)
  }
}
