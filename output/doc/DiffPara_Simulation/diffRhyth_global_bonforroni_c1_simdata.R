##differential rhythmicity: global test
#
#speed3
# Set directory and load libraries-----------------------------------------------------------
rm(list = ls())
dir = "/home/xix66/circadian/ThePipeline/diffRhyth_global_bonferroni_noR2change"
setwd(dir)

library(parallel)
#simulate data: each data will have the same set of the genes, just repeat the following for 10 times
for(r in 1:100){
  set.seed(r)
  n.genes = 3000
  n.samples.each = 20
  n.rhythmic.one = 1000
  n.rhythmic.both = 1000
  #parameters for only one rhythmic
  amp.one.1 = runif(n.rhythmic.one, min = 1, max = 3)
  phase.one.1 = runif(n.rhythmic.one, min = 0, max = 2*pi)
  offset.one.1 = runif(n.rhythmic.one, min = 5, max = 10)
  amp.one.2 = runif(n.rhythmic.one, min = 1, max = 3)
  phase.one.2 = runif(n.rhythmic.one, min = 0, max = 2*pi)
  offset.one.2 = runif(n.rhythmic.one, min = 5, max = 10)
  
  #parameters for both rhythmic 
  amp.both.2 = amp.both.1 = runif(n.rhythmic.both, min = 1, max = 3)
  phase.both.2 = phase.both.1 = runif(n.rhythmic.both, min = 0, max = 2*pi)
  offset.both.2 = offset.both.1 = runif(n.rhythmic.both, min = 5, max = 10)
  # sd.both.2 = sd.both.1 = rep(1, n.rhythmic.both)
  # sd.change.vec = c(0.5, 1.5)
  #change parameters
  ind.list = list(
    ind.nochange = 1:650,
    ind.change.offset = 651:700,
    ind.change.A = 701:750,
    ind.change.phase = 751:800,
    ind.change.A.phase = 801:850,
    ind.change.A.offset = 851:900,
    ind.change.phase.offset = 901:950,
    ind.change.A.phase.offset = 951:1000
  )
  #
  offset.both.2[do.call(c, ind.list[grepl("offset", names(ind.list))])] = 1.25*offset.both.2[do.call(c, ind.list[grepl("offset", names(ind.list))])] 
  amp.both.2[do.call(c, ind.list[grepl("A", names(ind.list))])] = 1.5*amp.both.2[do.call(c, ind.list[grepl("A", names(ind.list))])] 
  phase.both.2[do.call(c, ind.list[grepl("phase", names(ind.list))])] = phase.both.2[do.call(c, ind.list[grepl("phase", names(ind.list))])]  + 0.5*pi

  noise.mat = matrix(rnorm(n.samples.each*2*n.genes, sd = 1), ncol = n.samples.each*2)
  # noise.mat[1:n.rhythmic.both, (n.samples.each+1):(2*n.samples.each)] = 
  #   do.call(rbind, lapply(sd.both.2, function(a.sd){rnorm(n.samples.each, sd = a.sd)}))
    
  tod.1 = runif(n.samples.each, min = 0, max = 24)
  tod.2 = runif(n.samples.each, min = 0, max = 24) 
  
  noise.mat[1:n.rhythmic.both, 1:n.samples.each] = 
    noise.mat[1:n.rhythmic.both, 1:n.samples.each] + do.call(rbind, lapply(1:n.rhythmic.both, function(i){
      amp.both.1[i] * cos(2*pi/24*tod.1+ phase.both.1[i]) + offset.both.1[i]
    }))
  noise.mat[1:n.rhythmic.both, (n.samples.each+1):(2*n.samples.each)] = 
    noise.mat[1:n.rhythmic.both, (n.samples.each+1):(2*n.samples.each)] + do.call(rbind, lapply(1:n.rhythmic.both, function(i){
      amp.both.2[i] * cos(2*pi/24*tod.2+ phase.both.2[i]) + offset.both.2[i]
    }))
  noise.mat[(n.rhythmic.both+1):(n.rhythmic.both+n.rhythmic.one), 1:n.samples.each] = 
    noise.mat[(n.rhythmic.both+1):(n.rhythmic.both+n.rhythmic.one), 1:n.samples.each] +
    do.call(rbind, lapply(1:n.rhythmic.one, function(i){
      amp.one.1[i] * cos(2*pi/24*tod.1+ phase.one.1[i]) + offset.one.1[i]
    }))
  noise.mat[(n.rhythmic.both+n.rhythmic.one+1):(n.rhythmic.both+n.rhythmic.one*2), (n.samples.each+1):(2*n.samples.each)] = 
    noise.mat[(n.rhythmic.both+n.rhythmic.one+1):(n.rhythmic.both+n.rhythmic.one*2), (n.samples.each+1):(2*n.samples.each)] +
    do.call(rbind, lapply(1:n.rhythmic.both, function(i){
      amp.one.2[i] * cos(2*pi/24*tod.2+ phase.one.2[i]) + offset.one.2[i]
    }))
  
  save.image(paste0(dir, "/simdata/data_rep", r, ".Rdata"))
}

