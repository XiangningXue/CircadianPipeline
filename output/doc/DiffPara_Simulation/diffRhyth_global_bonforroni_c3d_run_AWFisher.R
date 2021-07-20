##differential rhythmicity: global test
#
#speed3
# Set directory and load libraries-----------------------------------------------------------
rm(list = ls())
dir = "/home/xix66/circadian/ThePipeline/diffRhyth_global_bonferroni_noR2change"
setwd(dir)


# Run with the bonferroni->BH->single test --------------------------------
library(parallel)
library(circacompare)
#simulate data: each data will have the same set of the genes, just repeat the following for 10 times
mclapply(1:100, function(r){
  circa.res.tab = readRDS(paste0(dir, "/res_circacompare/res", "_n", r, ".rds"))
  
  #do correction 
  AW.res = AWFisher_pvalue(as.matrix(circa.res.tab[, c("p.Mesor", "p.amp", "p.phase")]))
  
  circa.res.tab$p.BH = p.adjust(AW.res$pvalues, "BH")
  circa.res.tab$ind.global = circa.res.tab$p.BH <0.05
  circa.res.tab$ind.Mesor = circa.res.tab$p.BH <0.05&AW.res$weights[, 1] == 1
  circa.res.tab$ind.amp = circa.res.tab$p.BH <0.05&AW.res$weights[, 2] == 1
  circa.res.tab$ind.phase = circa.res.tab$p.BH <0.05&AW.res$weights[, 3] == 1

  out = circa.res.tab
  saveRDS(out, file = paste0(dir, "/res_circacompare_AWFisher/res", "_n", r, ".rds"))
}, mc.cores = 20)


