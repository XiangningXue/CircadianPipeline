diff_rhythmicity_bootstrap = function(x1 = x1, x2 = x2, subset = x.joint$Rhythmic.GT.One, period = period, B = nSampling,
                                        Sampling.save = Sampling.save, Sampling.file.label, parallel, cores){

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

  # x12.overlap = cbind(x1.overlap, x2.overlap)
  # t.all = c(x1$tod, x2$tod)
  t1 = x1$tod
  t2 = x2$tod

  #indexes = seq_len(length(t.all))
  index1.0 = seq_len(length(t1))
  index2.0 = seq_len(length(t2))
  # n1 = length(x1$tod)
  # n2 = length(x2$tod)

  boot.list = option_parallel(1:B, function(b){
    set.seed(b)
    index1 = sample(index1.0, replace = TRUE)
    index2 = sample(index2.0, replace = TRUE)

    one.boot.list = lapply(1:length(overlap.g), function(a){
      x1.boot.res = one_cosinor_OLS(t1[index1], y = as.numeric(x1.overlap[a, index1]), alpha, period, CI = FALSE)
      x1.boot.res.tab = as.data.frame(list(M = x1.boot.res$M$est,
                                           A = x1.boot.res$A$est,
                                           phase = x1.boot.res$phase$est, peak = x1.boot.res$peak,
                                           sigma2 = x1.boot.res$test$sigma2,
                                           R2 = x1.boot.res$test$R2))
      x2.boot.res = one_cosinor_OLS(t2[index2], y = as.numeric(x2.overlap[a, index2]), alpha, period, CI = FALSE)
      x2.boot.res.tab = as.data.frame(list(M = x2.boot.res$M$est,
                                           A = x2.boot.res$A$est,
                                           phase = x2.boot.res$phase$est, peak = x2.boot.res$peak,
                                           sigma2 = x2.boot.res$test$sigma2,
                                           R2 = x2.boot.res$test$R2))
      return(list(x1.boot = x1.boot.res.tab,
                  x2.boot = x2.boot.res.tab))
    })

    one.x1.boot.tab = do.call(rbind.data.frame, lapply(one.boot.list, `[[`, 1))
    one.x2.boot.tab = do.call(rbind.data.frame, lapply(one.boot.list, `[[`, 2))


    diffparas_null = list(M_null = one.x2.boot.tab$M-one.x1.boot.tab$M,
                          A_null = one.x2.boot.tab$A-one.x1.boot.tab$A,
                          phase_null = choose.abs.min(cbind(one.x2.boot.tab$phase-one.x1.boot.tab$phase,
                                                            one.x2.boot.tab$phase+2*pi-one.x1.boot.tab$phase,
                                                            one.x2.boot.tab$phase-2*pi-one.x1.boot.tab$phase)),
                          peak_null = choose.abs.min(cbind(one.x2.boot.tab$peak-one.x1.boot.tab$peak,
                                                           one.x2.boot.tab$peak+period-one.x1.boot.tab$peak,
                                                           one.x2.boot.tab$peak-period-one.x1.boot.tab$peak)),
                          R2_null = one.x2.boot.tab$R2-one.x1.boot.tab$R2)

    if(!is.null(Sampling.save)){
      save(one.x1.boot.tab, one.x2.boot.tab, index1, index2, overlap.g,
           file = paste0(file.path(Sampling.save, paste0(Sampling.file.label, "_Bootstrap", b, ".rData"))))
    }

    if(b%%(B/100)==0){
      cat(paste0("Differential Rhythmicity Analysis Bootstrap step: ", round(b/B, 2)*100, "% done. \n"))
    }

    return(diffparas_null)
  }, parallel, cores)

  M_null = do.call(cbind, lapply(boot.list, `[[`, "M_null")); M_null_sd = apply(M_null, 1, sd)
  A_null = do.call(cbind, lapply(boot.list, `[[`, "A_null")); A_null_sd = apply(A_null, 1, sd)
  phase_null = do.call(cbind, lapply(boot.list, `[[`, "phase_null")); phase_null_sd = apply(phase_null, 1, sd)
  peak_null = do.call(cbind, lapply(boot.list, `[[`, "peak_null")); peak_null_sd = apply(peak_null, 1, sd)
  R2_null = do.call(cbind, lapply(boot.list, `[[`, "R2_null")); R2_null_sd = apply(R2_null, 1, sd)

  x1.true = x1$rhythm[match(overlap.g, x1$label), ]
  x2.true = x2$rhythm[match(overlap.g, x2$label), ]

  M_obs = x2.true$M-x1.true$M
  A_obs = x2.true$A-x1.true$A
  phase_obs = choose.abs.min(cbind(x2.true$phase-x1.true$phase, x2.true$phase+2*pi-x1.true$phase, x2.true$phase-2*pi-x1.true$phase))
  peak_obs = choose.abs.min(cbind(x2.true$peak-x1.true$peak, x2.true$peak+period-x1.true$peak, x2.true$peak-period-x1.true$peak))
  R2_obs = x2.true$R2-x1.true$R2

  #ap_R2_perGene <- apply(R2_null - R2_obs,1,function(x) min(mean(x >= 0), mean(x <= 0)) * 2)

  ## Sampling p-value by pooling all genes together
  ap_M_allGene0 <- pnorm(M_obs/M_null_sd)
  ap_M_allGene <- pmin(ap_M_allGene0, 1 - ap_M_allGene0) * 2

  ap_A_allGene0 <- pnorm(A_obs/A_null_sd)
  ap_A_allGene <- pmin(ap_A_allGene0, 1 - ap_A_allGene0) * 2

  ap_phase_allGene0 <- pnorm(phase_obs/phase_null_sd)
  ap_phase_allGene <- pmin(ap_phase_allGene0, 1 - ap_phase_allGene0) * 2

  ap_peak_allGene0 <- pnorm(peak_obs/peak_null_sd)
  ap_peak_allGene <- pmin(ap_peak_allGene0, 1 - ap_peak_allGene0) * 2

  ap_R2_allGene0 <- pnorm(R2_obs/R2_null_sd)
  ap_R2_allGene <- pmin(ap_R2_allGene0, 1 - ap_R2_allGene0) * 2

  out = data.frame(label = overlap.g,
                   delta.M = M_obs,
                   p.M = ap_M_allGene,
                   delta.A = A_obs,
                   p.A = ap_A_allGene,
                   delta.phase = phase_obs,
                   delta.peak = peak_obs,
                   p.phase = ap_phase_allGene,
                   #p.peak = ap_peak_allGene, #p for bootutaiton p for peak and phase are almost the same, no need for both
                   delta.R2 = R2_obs,
                   p.R2 = ap_R2_allGene)
  if(!is.null(Sampling.save)){
    write.csv(out, paste0(file.path(Sampling.save, paste0(Sampling.file.label, "_Bootstrap", "_out", ".csv"))))
  }
  return(out)
}
