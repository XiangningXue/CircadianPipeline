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
  p.vec.allthree = c(circa.res.tab$p.Mesor, circa.res.tab$p.amp, circa.res.tab$p.phase)
  q.vec.allthree = p.adjust(p.vec.allthree, "BH")
  circa.res.tab$q.Mesor = q.vec.allthree[1:nrow(circa.res.tab)]
  circa.res.tab$q.amp = q.vec.allthree[(nrow(circa.res.tab)+1):(nrow(circa.res.tab)*2)]
  circa.res.tab$q.phase = q.vec.allthree[(nrow(circa.res.tab)*2+1):(nrow(circa.res.tab)*3)]
  circa.res.tab$ind.global = 3*pmin(circa.res.tab$q.Mesor, circa.res.tab$q.amp, circa.res.tab$q.phase)<0.05
  circa.res.tab$ind.Mesor = circa.res.tab$q.Mesor < 0.05/3
  circa.res.tab$ind.amp = circa.res.tab$q.amp < 0.05/3
  circa.res.tab$ind.phase = circa.res.tab$q.phase < 0.05/3

  out = circa.res.tab
  saveRDS(out, file = paste0(dir, "/res_circacompare_BH/res", "_n", r, ".rds"))
}, mc.cores = 20)


