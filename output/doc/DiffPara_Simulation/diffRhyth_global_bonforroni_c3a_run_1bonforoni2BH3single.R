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
n.rep = seq(1, 10)
mclapply(1:100, function(r){
  circa.res.tab = readRDS(paste0(dir, "/res_circacompare/res", "_n", r, ".rds"))
  #do correction 
  circa.res.tab$p.InGeneCorrected = pmin(3*pmin(circa.res.tab$p.Mesor, circa.res.tab$p.amp, circa.res.tab$p.phase), 1)
  circa.res.tab$p.BH = p.adjust(circa.res.tab$p.InGeneCorrected, "BH")
  the.BH.p.cutoff = max(circa.res.tab$p.InGeneCorrected[circa.res.tab$p.BH<0.05])
  circa.res.tab$ind.global = circa.res.tab$p.BH <0.05
  circa.res.tab$ind.Mesor = circa.res.tab$p.Mesor<the.BH.p.cutoff/3
  circa.res.tab$ind.amp = circa.res.tab$p.amp<the.BH.p.cutoff/3
  circa.res.tab$ind.phase = circa.res.tab$p.phase<the.BH.p.cutoff/3
  
  out = circa.res.tab
  saveRDS(out, file = paste0(dir, "/res_circacompare_1Bonferroni2BH3Individual/res", "_n", r, ".rds"))
}, mc.cores = 20)


