##differential rhythmicity: global test
#
#speed3
# Set directory and load libraries-----------------------------------------------------------
rm(list = ls())
dir = "/home/xix66/circadian/ThePipeline/diffRhyth_global_bonferroni_noR2change"
setwd(dir)


# Run with the bonferroni->BH->single test --------------------------------
library(parallel)
#simulate data: each data will have the same set of the genes, just repeat the following for 10 times
system(paste0("mkdir -p ", paste0(dir, "/res_cosinor_permutation/res_permutation_1Bonferroni2BH3Individual")))
res1 = mclapply(1:100, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  #do correction 
  rthmicChange$p.InGeneCorrected = pmin(3*pmin(rthmicChange$intshiftPvalue, rthmicChange$AshiftPvalue, rthmicChange$phasePvalue), 1)
  rthmicChange$p.BH = p.adjust(rthmicChange$p.InGeneCorrected, "BH")
  the.BH.p.cutoff = max(rthmicChange$p.InGeneCorrected[rthmicChange$p.BH<0.05])
  rthmicChange$ind.global = rthmicChange$p.BH <0.05
  rthmicChange$ind.Mesor = rthmicChange$intshiftPvalue<the.BH.p.cutoff/3
  rthmicChange$ind.amp = rthmicChange$AshiftPvalue<the.BH.p.cutoff/3
  rthmicChange$ind.phase = rthmicChange$phasePvalue<the.BH.p.cutoff/3
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_1Bonferroni2BH3Individual/res", "_n", r, ".rds"))
}, mc.cores = 10)
which(sapply(res1, function(a){length(a)})==1)
res1.2 = mclapply(38, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  #do correction 
  rthmicChange$p.InGeneCorrected = pmin(3*pmin(rthmicChange$intshiftPvalue, rthmicChange$AshiftPvalue, rthmicChange$phasePvalue), 1)
  rthmicChange$p.BH = p.adjust(rthmicChange$p.InGeneCorrected, "BH")
  the.BH.p.cutoff = max(rthmicChange$p.InGeneCorrected[rthmicChange$p.BH<0.05])
  rthmicChange$ind.global = rthmicChange$p.BH <0.05
  rthmicChange$ind.Mesor = rthmicChange$intshiftPvalue<the.BH.p.cutoff/3
  rthmicChange$ind.amp = rthmicChange$AshiftPvalue<the.BH.p.cutoff/3
  rthmicChange$ind.phase = rthmicChange$phasePvalue<the.BH.p.cutoff/3
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_1Bonferroni2BH3Individual/res", "_n", r, ".rds"))
}, mc.cores = 10)
which(sapply(res1.2, function(a){length(a)})==1)

# Run with the bonferroni->BH->single test --------------------------------
system(paste0("mkdir -p ", paste0(dir, "/res_cosinor_permutation/res_permutation_1BH2Bonferroni")))
res2 = mclapply(1:100, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  
  #do correction 
  rthmicChange$q.Mesor = p.adjust(rthmicChange$intshiftPvalue, "BH")
  rthmicChange$q.amp = p.adjust(rthmicChange$AshiftPvalue, "BH")
  rthmicChange$q.phase = p.adjust(rthmicChange$phasePvalue, "BH")
  rthmicChange$ind.global = 3*pmin(rthmicChange$q.Mesor, rthmicChange$q.amp, rthmicChange$q.phase)<0.05
  rthmicChange$ind.Mesor = rthmicChange$q.Mesor < 0.05/3
  rthmicChange$ind.amp = rthmicChange$q.amp < 0.05/3
  rthmicChange$ind.phase = rthmicChange$q.phase < 0.05/3
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_1BH2Bonferroni/res", "_n", r, ".rds"))
}, mc.cores = 20)

# Run with the bonferroni->BH->single test --------------------------------
system(paste0("mkdir -p ", paste0(dir, "/res_cosinor_permutation/res_permutation_1BH2Bonferroni0")))
res3 = mclapply(1:100, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  
  #do correction 
  rthmicChange$q.Mesor = p.adjust(rthmicChange$intshiftPvalue, "BH")
  rthmicChange$q.amp = p.adjust(rthmicChange$AshiftPvalue, "BH")
  rthmicChange$q.phase = p.adjust(rthmicChange$phasePvalue, "BH")
  rthmicChange$ind.global = 3*pmin(rthmicChange$q.Mesor, rthmicChange$q.amp, rthmicChange$q.phase)<0.05
  rthmicChange$ind.Mesor = rthmicChange$q.Mesor < 0.05
  rthmicChange$ind.amp = rthmicChange$q.amp < 0.05
  rthmicChange$ind.phase = rthmicChange$q.phase < 0.05
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_1BH2Bonferroni0/res", "_n", r, ".rds"))
}, mc.cores = 20)

# Run with the bonferroni->BH->single test --------------------------------
system(paste0("mkdir -p ", paste0(dir, "/res_cosinor_permutation/res_permutation_BH")))
res4 = mclapply(1:100, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  
  #do correction 
  p.vec.allthree = c(rthmicChange$intshiftPvalue, rthmicChange$AshiftPvalue, rthmicChange$phasePvalue)
  q.vec.allthree = p.adjust(p.vec.allthree, "BH")
  rthmicChange$q.Mesor = q.vec.allthree[1:nrow(rthmicChange)]
  rthmicChange$q.amp = q.vec.allthree[(nrow(rthmicChange)+1):(nrow(rthmicChange)*2)]
  rthmicChange$q.phase = q.vec.allthree[(nrow(rthmicChange)*2+1):(nrow(rthmicChange)*3)]
  rthmicChange$ind.global = 3*pmin(rthmicChange$q.Mesor, rthmicChange$q.amp, rthmicChange$q.phase)<0.05
  rthmicChange$ind.Mesor = rthmicChange$q.Mesor < 0.05/3
  rthmicChange$ind.amp = rthmicChange$q.amp < 0.05/3
  rthmicChange$ind.phase = rthmicChange$q.phase < 0.05/3
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_BH/res", "_n", r, ".rds"))
}, mc.cores = 20)

# Run with the bonferroni->BH->single test --------------------------------
system(paste0("mkdir -p ", paste0(dir, "/res_cosinor_permutation/res_permutation_BH0")))
res5 = mclapply(1:100, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  
  #do correction 
  p.vec.allthree = c(rthmicChange$intshiftPvalue, rthmicChange$AshiftPvalue, rthmicChange$phasePvalue)
  q.vec.allthree = p.adjust(p.vec.allthree, "BH")
  rthmicChange$q.Mesor = q.vec.allthree[1:nrow(rthmicChange)]
  rthmicChange$q.amp = q.vec.allthree[(nrow(rthmicChange)+1):(nrow(rthmicChange)*2)]
  rthmicChange$q.phase = q.vec.allthree[(nrow(rthmicChange)*2+1):(nrow(rthmicChange)*3)]
  rthmicChange$ind.global = 3*pmin(rthmicChange$q.Mesor, rthmicChange$q.amp, rthmicChange$q.phase)<0.05
  rthmicChange$ind.Mesor = rthmicChange$q.Mesor < 0.05
  rthmicChange$ind.amp = rthmicChange$q.amp < 0.05
  rthmicChange$ind.phase = rthmicChange$q.phase < 0.05
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_BH0/res", "_n", r, ".rds"))
}, mc.cores = 20)

# Run with the bonferroni->BH->single test --------------------------------
system(paste0("mkdir -p ", paste0(dir, "/res_cosinor_permutation/res_permutation_AWFisher")))
library(AWFisher)
res6 = mclapply(1:100, function(r){
  rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
  
  #do correction 
  AW.res = AWFisher_pvalue(as.matrix(rthmicChange[, c("intshiftPvalue", "AshiftPvalue", "phasePvalue")]))
  
  rthmicChange$p.BH = p.adjust(AW.res$pvalues, "BH")
  rthmicChange$ind.global = rthmicChange$p.BH <0.05
  rthmicChange$ind.Mesor = rthmicChange$p.BH <0.05&AW.res$weights[, 1] == 1
  rthmicChange$ind.amp = rthmicChange$p.BH <0.05&AW.res$weights[, 2] == 1
  rthmicChange$ind.phase = rthmicChange$p.BH <0.05&AW.res$weights[, 3] == 1
  
  out = rthmicChange
  saveRDS(out, file = paste0(dir, "/res_cosinor_permutation/res_permutation_AWFisher/res", "_n", r, ".rds"))
}, mc.cores = 20)

