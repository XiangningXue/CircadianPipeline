rhythmicity_permutation = function(X = X, tod = tod, period = period, B = nPermutation,
                                   permutation.save = permutation.save, permutation.file.label = permutation.file.label, parallel, cores){
  n = nrow(X)

  true.list = option_parallel(1:nrow(X), function(a){
    true.res = one_cosinor_OLS(tod, y = as.numeric(X[a, ]), alpha, period, CI = FALSE)
    true.res.tab = as.data.frame(list(M = true.res$M$est,
                                      A = true.res$A$est,
                                      phase = true.res$phase$est, peak = true.res$peak,
                                      sigma2 = true.res$test$sigma2,
                                      R2 = true.res$test$R2))

    return(true.res.tab)
  }, parallel, cores)

  perm.list = option_parallel(1:B, function(b){
    set.seed(b)
    shuffleTOD = sample(tod)
    one.perm.list = lapply(1:nrow(X), function(a){
      perm.res = one_cosinor_OLS(shuffleTOD, y = as.numeric(X[a, ]), alpha, period, CI = FALSE)
      perm.res.tab = as.data.frame(list(M = perm.res$M$est,
                                        A = perm.res$A$est,
                                        phase = perm.res$phase$est, peak = perm.res$peak,
                                        sigma2 = perm.res$test$sigma2,
                                        R2 = perm.res$test$R2))

      return(perm.res.tab)
    })
    one.perm.tab = do.call(rbind.data.frame, one.perm.list)

    if(permutation.save!="NULL"){
      saveRDS(one.perm.tab, paste0(file.path(permutation.save, paste0(permutation.file.label, "_PermWithin", b, ".rds"))))
    }

    if(b%%(B/100)==0){
      cat(paste0("Rhythmicity Analysis permutation step: ", round(b/B, 2)*100, "% done. \n"))
    }
    return(one.perm.tab)
  }, parallel, cores)

  res.true.tab = do.call(rbind.data.frame, true.list)
  R2.pool = c(res.true.tab$R2, unlist(lapply(perm.list, function(a){a$R2})))
  R2.Rank <- 1 - (rank(R2.pool)[1:length(res.true.tab$R2)] - 0.5)/length(R2.pool) #the 0.5 is to avoid p-value = 0
  res.true.tab$pvalue <- R2.Rank


  return(res.true.tab)
}


