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
for(r in 1:100){
  load(paste0(dir, "/simdata/data_rep", r, ".Rdata"))
  #now we will only focus on the first 10000 genes. #both rhythmic
  dat = noise.mat[1:1000, ]
  
  #run circacompare
  exp.list = apply(dat, 1, function(one.gene){
    one.gene.data = data.frame(time = c(tod.1, tod.2), measure = one.gene, group = factor(c(rep(1, length(tod.1)), rep(2, length(tod.2)))))
    return(one.gene.data)
  })
  circa.res = mclapply(exp.list, function(one.data){
    one.res = circacompare(one.data, col_time = "time", col_group = "group", col_outcome = "measure", 
                           alpha_threshold = 1)
    return(as.data.frame(list(rhyth.1 = one.res[[2]][2, 2], 
                              rhyth.2 = one.res[[2]][3, 2],
                              p.Mesor = one.res[[2]][7, 2],
                              p.amp = one.res[[2]][11, 2],
                              p.phase = one.res[[2]][15, 2])))
  }, mc.cores = 20)
  
  circa.res.tab = do.call(rbind.data.frame, circa.res)
  
  out = circa.res.tab
  saveRDS(out, file = paste0(dir, "/res_circacompare/res", "_n", r, ".rds"))
}


