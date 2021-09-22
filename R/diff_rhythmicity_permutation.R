diff_rhythmicity_permutation = function(x1 = x1, x2 = x2, subset = x.joint$Rhythmic.GT.One, period = period, B = nPermutation,
                                        permutation.save = permutation.save, permutation.file.label, parallel, cores){

  choose.abs.min = function(mat = matrix(rnorm(21), ncol = 3)){
    abs.min.vec = apply(mat, 1, function(a){
      min.abs.idx = which.min(abs(a))
      return(a[min.abs.idx])
    })
    return(abs.min.vec)
  }

  overlap.g = subset

  x1.overlap = x1$data[match(overlap.g, x1$label), ]
  x2.overlap = x2$data[match(overlap.g, x2$label), ]

  x12.overlap = cbind(x1.overlap, x2.overlap)
  t.all = c(x1$tod, x2$tod)

  indexes = seq_len(length(t.all))
  n1 = length(x1$tod)

  perm.list = option_parallel(1:B, function(b){
    set.seed(b)
    index1 = sample(indexes, n1)
    index2 = setdiff(indexes, index1)

    one.perm.list = lapply(1:nrow(x12.overlap), function(a){
      x1.perm.res = one_cosinor_OLS(t.all[index1], y = as.numeric(x12.overlap[a, index1]), alpha, period, CI = FALSE)
      x1.perm.res.tab = as.data.frame(list(M = x1.perm.res$M$est,
                                           A = x1.perm.res$A$est,
                                           phase = x1.perm.res$phase$est, peak = x1.perm.res$peak,
                                           sigma2 = x1.perm.res$test$sigma2,
                                           R2 = x1.perm.res$test$R2))
      x2.perm.res = one_cosinor_OLS(t.all[index2], y = as.numeric(x12.overlap[a, index2]), alpha, period, CI = FALSE)
      x2.perm.res.tab = as.data.frame(list(M = x2.perm.res$M$est,
                                           A = x2.perm.res$A$est,
                                           phase = x2.perm.res$phase$est, peak = x2.perm.res$peak,
                                           sigma2 = x2.perm.res$test$sigma2,
                                           R2 = x2.perm.res$test$R2))
      return(list(x1.perm = x1.perm.res.tab,
                  x2.perm = x2.perm.res.tab))
    })

    one.x1.perm.tab = do.call(rbind.data.frame, lapply(one.perm.list, `[[`, 1))
    one.x2.perm.tab = do.call(rbind.data.frame, lapply(one.perm.list, `[[`, 2))


    diffparas_null = list(M_null = one.x2.perm.tab$M-one.x1.perm.tab$M,
                          A_null = one.x2.perm.tab$A-one.x1.perm.tab$A,
                          phase_null = choose.abs.min(cbind(one.x2.perm.tab$phase-one.x1.perm.tab$phase,
                                                            one.x2.perm.tab$phase+2*pi-one.x1.perm.tab$phase,
                                                            one.x2.perm.tab$phase-2*pi-one.x1.perm.tab$phase)),
                          peak_null = choose.abs.min(cbind(one.x2.perm.tab$peak-one.x1.perm.tab$peak,
                                                           one.x2.perm.tab$peak+period-one.x1.perm.tab$peak,
                                                           one.x2.perm.tab$peak-period-one.x1.perm.tab$peak)),
                          R2_null = one.x2.perm.tab$R2-one.x1.perm.tab$R2)

    if(!is.null(permutation.save)){
      save(one.x1.perm.tab, one.x2.perm.tab, index1, index2, overlap.g,
           file = paste0(file.path(permutation.save, paste0(permutation.file.label, "_PermBetween", b, ".rData"))))
    }

    if(b%%(B/100)==0){
      cat(paste0("Differential Rhythmicity Analysis permutation step: ", round(b/B, 2)*100, "% done. \n"))
    }

    return(diffparas_null)
  }, parallel, cores)

  M_null = do.call(cbind, lapply(perm.list, `[[`, "M_null"))
  A_null = do.call(cbind, lapply(perm.list, `[[`, "A_null"))
  phase_null = do.call(cbind, lapply(perm.list, `[[`, "phase_null"))
  peak_null = do.call(cbind, lapply(perm.list, `[[`, "peak_null"))
  R2_null = do.call(cbind, lapply(perm.list, `[[`, "R2_null"))

  x1.true = x1$rhythm[match(overlap.g, x1$label), ]
  x2.true = x2$rhythm[match(overlap.g, x2$label), ]

  M_obs = x2.true$M-x1.true$M
  A_obs = x2.true$A-x1.true$A
  phase_obs = choose.abs.min(cbind(x2.true$phase-x1.true$phase, x2.true$phase+2*pi-x1.true$phase, x2.true$phase-2*pi-x1.true$phase))
  peak_obs = choose.abs.min(cbind(x2.true$peak-x1.true$peak, x2.true$peak+period-x1.true$peak, x2.true$peak-period-x1.true$peak))
  R2_obs = x2.true$R2-x1.true$R2

  #ap_R2_perGene <- apply(R2_null - R2_obs,1,function(x) min(mean(x >= 0), mean(x <= 0)) * 2)

  ## permutation p-value by pooling all genes together
  ap_M_allGene0 <- rank(cbind(M_obs, M_null))[1:length(overlap.g)]/(length(overlap.g)*(B + 1))
  ap_M_allGene <- pmin(ap_M_allGene0, 1 - ap_M_allGene0) * 2

  ap_A_allGene0 <- rank(cbind(A_obs, A_null))[1:length(overlap.g)]/(length(overlap.g)*(B + 1))
  ap_A_allGene <- pmin(ap_A_allGene0, 1 - ap_A_allGene0) * 2

  ap_phase_allGene0 <- rank(cbind(phase_obs, phase_null))[1:length(overlap.g)]/(length(overlap.g)*(B + 1))
  ap_phase_allGene <- pmin(ap_phase_allGene0, 1 - ap_phase_allGene0) * 2

  ap_peak_allGene0 <- rank(cbind(peak_obs, peak_null))[1:length(overlap.g)]/(length(overlap.g)*(B + 1))
  ap_peak_allGene <- pmin(ap_peak_allGene0, 1 - ap_peak_allGene0) * 2

  ap_R2_allGene0 <- rank(cbind(R2_obs, R2_null))[1:length(overlap.g)]/(length(overlap.g)*(B + 1))
  ap_R2_allGene <- pmin(ap_R2_allGene0, 1 - ap_R2_allGene0) * 2

  out = data.frame(label = overlap.g,
                   delta.M = M_obs,
                   p.M = ap_M_allGene,
                   delta.A = A_obs,
                   p.A = ap_A_allGene,
                   delta.phase = phase_obs,
                   delta.peak = peak_obs,
                   p.phase = ap_phase_allGene,
                   #p.peak = ap_peak_allGene, #p for permutaiton p for peak and phase are almost the same, no need for both
                   delta.R2 = R2_obs,
                   p.R2 = ap_R2_allGene)
  if(!is.null(permutation.save)){
    write.csv(out, paste0(file.path(permutation.save, paste0(permutation.file.label, "_PermBetween", "_out", ".csv"))))
  }
  return(out)
}
