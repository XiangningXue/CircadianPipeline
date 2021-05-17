option_parallel = function(X=x, FUN=fun, parallel = TRUE, cores = 1){
  if(parallel){
    parallel::mclapply(X, FUN, mc.cores = cores)
  }else{
    lapply(X, FUN)
  }
}
