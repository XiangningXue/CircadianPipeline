adj_minP = function(diffPar.tab = diffPar.tab, p.adjust.method = "BH", alpha = 0.05){

  diffPar.tab$p.combined = pmin(3*pmin(diffPar.tab$p.M, diffPar.tab$p.A, diffPar.tab$p.phase), 1)
  diffPar.tab$q.combined = stats::p.adjust(diffPar.tab$p.combined, p.adjust.method)
  the.BH.p.cutoff = max(diffPar.tab$p.combined[diffPar.tab$q.combined<alpha])
  diffPar.tab$ind.global = diffPar.tab$q.combined <alpha
  diffPar.tab$ind.M = diffPar.tab$p.M<the.BH.p.cutoff/3
  diffPar.tab$ind.A = diffPar.tab$p.A<the.BH.p.cutoff/3
  diffPar.tab$ind.phase = diffPar.tab$p.phase<the.BH.p.cutoff/3

  return(diffPar.tab)
}


adj_BH_separate = function(diffPar.tab = diffPar.tab, p.adjust.method = "BH", alpha = 0.05){

  diffPar.tab$q.M = stats::p.adjust(diffPar.tab$p.M, p.adjust.method)
  diffPar.tab$q.A = stats::p.adjust(diffPar.tab$p.A, p.adjust.method)
  diffPar.tab$q.phase = stats::p.adjust(diffPar.tab$p.phase, p.adjust.method)
  diffPar.tab$ind.global = pmin(diffPar.tab$q.M, diffPar.tab$q.A, diffPar.tab$q.phase)<alpha
  diffPar.tab$ind.M = diffPar.tab$q.M<alpha
  diffPar.tab$ind.A = diffPar.tab$q.A<alpha
  diffPar.tab$ind.phase = diffPar.tab$q.phase<alpha

  return(diffPar.tab)
}

adj_BH_pool = function(diffPar.tab = diffPar.tab, p.adjust.method = "BH", alpha = 0.05){

  p.pool = c(diffPar.tab$p.M, diffPar.tab$p.A, diffPar.tab$p.phase)
  q.pool = stats::p.adjust(p.pool)
  diffPar.tab$q.M = q.pool[1:nrow(diffPar.tab)]
  diffPar.tab$q.A = q.pool[(nrow(diffPar.tab)+1):(nrow(diffPar.tab)*2)]
  diffPar.tab$q.phase = q.pool[(2*nrow(diffPar.tab)+1):(nrow(diffPar.tab)*3)]
  diffPar.tab$ind.global = pmin(diffPar.tab$q.M, diffPar.tab$q.A, diffPar.tab$q.phase)<alpha
  diffPar.tab$ind.M = diffPar.tab$q.M<alpha
  diffPar.tab$ind.A = diffPar.tab$q.A<alpha
  diffPar.tab$ind.phase = diffPar.tab$q.phase<alpha

  return(diffPar.tab)
}

adj_AWFisher = function(diffPar.tab = diffPar.tab, p.adjust.method = "BH", alpha = 0.05){
  #do correction
  AW.res = AWFisher::AWFisher_pvalue(as.matrix(diffPar.tab[, c("p.M", "p.A", "p.phase")]))

  diffPar.tab$p.combined = AW.res$pvalues
  diffPar.tab$q.combined = stats::p.adjust(AW.res$pvalues, p.adjust.method)
  diffPar.tab$ind.global = diffPar.tab$q.combined < alpha
  diffPar.tab$ind.M = diffPar.tab$q.combined <0.05&AW.res$weights[, 1] == 1
  diffPar.tab$ind.A = diffPar.tab$q.combined <0.05&AW.res$weights[, 2] == 1
  diffPar.tab$ind.phase = diffPar.tab$q.combined <0.05&AW.res$weights[, 3] == 1
  return(diffPar.tab)
}
