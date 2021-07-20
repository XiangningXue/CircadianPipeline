##differential rhythmicity: global test
#
#speed3
# Set directory and load libraries-----------------------------------------------------------
rm(list = ls())
dir = "/home/xix66/circadian/ThePipeline/diffRhyth_global_bonferroni_noR2change"
source("/home/xix66/circadian/ThePipeline/diffR2/cosinor package_CI2.R")
setwd(dir)


# Run with the bonferroni->BH->single test --------------------------------
library(parallel)
library(circacompare)
#simulate data: each data will have the same set of the genes, just repeat the following for 10 times
n.rep = seq(1, 10)
for(r in 1:100){
  
  
  # res.exist = list.files(paste0(dir, "/res_cosinor_permutation"))
  # if(sum(grepl(paste0("res", "_n", r, ".rds"), res.exist))){
  #   next
  # }
  
  load(paste0(dir, "/simdata/data_rep", r, ".Rdata"))
  #now we will only focus on the first 10000 genes. #both rhythmic
  dat = noise.mat[1:1000, ]
  
  #run cosinor 
  groups = factor(c(rep(1, length(tod.1)), rep(2, length(tod.2))))
  x.cov.levels = unique(groups)
  tod = c(tod.1, tod.2)
  rhyth.res = lapply(1:length(x.cov.levels), function(l){
    one.x.col.level = x.cov.levels[l]
    dat.sub = dat[, groups==one.x.col.level]
    print(dim(dat.sub))
    tod.sub = tod[groups==one.x.col.level]
    out.list = lapply(1:nrow(dat.sub), function(i){
      res = one_cosinor_OLS(tod = tod.sub, y = dat.sub[i, ])
      res.onerow = as.data.frame(list(offset = res$M$est,
                                      A = res$A$est,
                                      phase = res$phase$est, 
                                      pvalue = res$test$pval, R2 = res$test$R2))
      return(res.onerow)
    })
    out.tab = do.call(rbind.data.frame, out.list)
    return(out.tab)
  })
  
  null.dir = paste0(dir, "/res_cosinor_permutation/Null_data")

  B = 200
  indexes = seq_len(ncol(dat))
  rhyth.perm = mclapply(1:B, function(b){
    set.seed(b)
    print(b)
    
    index1 <- sample(indexes,length(tod.1))
    index2 <- setdiff(indexes, index1)  
    
    rhyth.perm1 = lapply(list(index1, index2), function(ind){
      out.list = lapply(1:nrow(dat), function(i){
        res = one_cosinor_OLS(tod = tod, y = dat[i,ind])
        res.onerow = as.data.frame(list(offset = res$M$est,
                                        A = res$A$est,
                                        phase = res$phase$est, 
                                        pvalue = res$test$pval, R2 = res$test$R2))
        return(res.onerow)
      })
      out.tab = do.call(rbind.data.frame, out.list)
      return(out.tab)
    })
    return(rhyth.perm1)
  }, mc.cores = 20)
  saveRDS(rhyth.perm, paste0(null.dir, "/perm_res_", r, ".rds"))
  
  #the true diff
  rhy.compare = rhyth.res[[2]]
  rhy.base = rhyth.res[[1]]
  R2change = rhy.compare$R2 - rhy.base$R2 #if is greater, than gain
  gain.sig = rhy.compare$pvalue<0.05 #if want to see significant gain, the case rhyth should be significant
  loss.sig = rhy.base$pvalue<0.05
  Ashift <- abs(abs(rhy.compare$A) - abs(rhy.base$A))
  phaseshift0 <- rhy.compare$phase - rhy.base$phase
  phaseshift <- pmin(abs(phaseshift0), 2*pi - abs(phaseshift0))
  intshift <- abs(rhy.compare$offset - rhy.base$offset)
  
  R2changeNULL <- R2change
  AshiftNULL <- Ashift
  phaseshiftNULL <- phaseshift
  intshiftNULL <- intshift
  
  #the permutation results
  res = mclapply(1:B, function(b){
    if(b%%50 == 0) print(b)
    null_para_base <- rhyth.perm[[b]][[1]]
    null_para_compare <- rhyth.perm[[b]][[2]]
    
    aR2change <- null_para_compare$R2 - null_para_base$R2
    aAshift <- abs(abs(null_para_compare$A) - abs(null_para_base$A))
    aphaseshift <- pmin(abs(null_para_compare$phase - null_para_base$phase),2*pi - abs(null_para_compare$phase - null_para_base$phase))
    aintshift <- abs(null_para_compare$offset - null_para_base$offset)
    
    return(list(aR2change  = aR2change, 
                aAshift = aAshift, 
                aphaseshift = aphaseshift, 
                aintshift = aintshift))
  }, mc.cores = 20)
  
  p <- nrow(rhy.base)
  
  R2changeNULL = append(R2change, do.call(c, lapply(res, function(b){b$aR2change})))
  AshiftNULL = append(Ashift, do.call(c, lapply(res, function(b){b$aAshift})))
  phaseshiftNULL = append(phaseshift, do.call(c, lapply(res, function(b){b$aphaseshift})))
  intshiftNULL = append(intshift, do.call(c, lapply(res, function(b){b$aintshift})))
  
  R2gainPvalue <- 1-rank(R2changeNULL)[1:p]/length(R2changeNULL) ## R2 female > male
  R2losePvalue <- rank(R2changeNULL)[1:p]/length(R2changeNULL) ## R2 female < male
  AshiftPvalue <- 1 - rank(AshiftNULL)[1:p]/length(AshiftNULL)
  phasePvalue <- 1 - rank(phaseshiftNULL)[1:p]/length(phaseshiftNULL)
  intshiftPvalue <- 1 - rank(intshiftNULL)[1:p]/length(intshiftNULL)
  R2changePvalue0 =  pmin(R2gainPvalue, 1 - R2losePvalue) * 2
  round(quantile(R2changeNULL[-c(1:p)]), 2)
  round(quantile(R2change), 2)
  
  R1 = do.call(c, lapply(rhyth.perm, function(b){b[[1]]$R2}))
  round(quantile(R1), 2)
  R2 = do.call(c, lapply(rhyth.perm, function(b){b[[2]]$R2}))
  round(quantile(R2), 2)
  
  
  sum(R2losePvalue<0.05)
  sum(R2losePvalue<0.05)
  sum(R2gainPvalue<0.05) 
  
  rthmicChange <- data.frame(genes=rownames(rhy.base),
                             rhyth.1 = rhy.base$pvalue, 
                             rhyth.2 = rhy.compare$pvalue,
                             R2.1 = rhy.base$R2,
                             R2.2 = rhy.compare$R2,
                             R2change = R2change, 
                             R2losePvalue=R2losePvalue,
                             R2gainPvalue=R2gainPvalue,
                             R2changePvalue0 = R2changePvalue0,
                             Ashift = Ashift, 
                             AshiftPvalue=AshiftPvalue,
                             phaseshift = phaseshift, 
                             phasePvalue=phasePvalue,
                             intshift = intshift,
                             intshiftPvalue=intshiftPvalue, 
                             gainSig = gain.sig, 
                             lossSig = loss.sig)
  
  saveRDS(rthmicChange, file = paste0(dir, "/res_cosinor_permutation/res", "_n", r, ".rds"))
}


