#' Title Rhythmicity Analysis for one data set
#'
#'This function fit the expression and time to this function \deqn{y = M+A*cos(2\pi/24*t+phase)+\epsilon}, where \eqn{\epsilon \sim N(0, \sigma^2)}.
#'
#' @param x  A list containing the following component: data: dataframe with rows of genes and columns of samples; tod: time of death or time of expression corresponding to the columns in data component; label: the gene names or other labels of the genes
#' @param period The length of the rhythmicity cycle. When it is 24, the signal is circadian.
#' @param method choose "OLS" or "permutation". OLS assumes Gaussian residual while permutation does not. Only OLS methods return confidence interval estimates. Notice that if permutation is selected, it is recommanded to use parallel computing.
#' @param nPermutation A numeric value. Number of permutation performer, only required when method = "permutation". The smallest possible p-value in the result will be 1/(n_gene*nPermutation), but larger nPermutation takes longer computing time.
#' @param permutation.save A local directory where you want to save the permutation result for future use. If "NULL" then no result will be saved.
#' @param permutation.file.label Special label for file name of permutation if you want to save the permutation file
#' @param alpha a number between 0 to 1. The critical level for the confidence interval
#' @param p.adjust.method choose valid input for p.adjust.method by checking p.adjust.methods for p.adjust() in R package `stat`
#' @param parallel TRUE/FALSE. If TRUE, parallel computing using mclapply will be used, which does not work on windows system.
#' @param cores number of cores used if using parallel computing
#'
#' @return A rhythm object, which is a list of original components in input x and rhythmicity results
#' @export
#'
#' @examples
#' #Simulate a time series data with 20 time points and 100 genes
#'
#' x.time = runif(20, min = 0, max = 24)
#' m = rnorm(100, 5); A = rnorm(100, 3); phase = runif(100, min = 0, max = 2*pi); sigma = 1
#' noise.mat = matrix(rnorm(100*20, 0, sigma), ncol = 20, nrow = 100)
#' signal.mat = t(sapply(1:100, function(a){m[a]+A[a]*cos(2*pi/24*x.time+phase[a])}))
#' x = list(data = as.data.frame(noise.mat + signal.mat),
#' tod = x.time,
#' label = paste("gene", seq_len(100)))
#' x.Rhythm = CP_Rhythmicity(x, parallel = FALSE)
#'
#'
CP_Rhythmicity = function(x = list(data = data1, tod = tod, label = label),
                          period = 24, method = "OLS", nPermutation=1000, permutation.save = getwd(), permutation.file.label = "Group1",
                          alpha = 0.05, p.adjust.method = "BH", parallel = TRUE, cores = 5){
  X = x$data
  t = x$tod
  g = x$label
  if(!is.data.frame(X)){
    stop("data must be a dataframe")
  }
  if(length(X)!=length(t)){
    stop("Number of samples in data does not match that in tod. ")
  }
  if(nrow(X)!=length(g)){
    stop("Number of labels does not match number of feature in data. ")
  }
  if(!is.numeric(t)){
    stop("Time must be numeric")
  }

  if(method == "OLS"){
    rhythmicity.list = option_parallel(1:nrow(X), function(a){
      res = one_cosinor_OLS(tod = t, y = as.numeric(X[a, ]), alpha, period, CI = TRUE)
      res.onerow = as.data.frame(list(label = g[a], P = period,
                                      M = res$M$est, A = res$A$est, phase = res$phase$est, peak = res$peak, pvalue = res$test$pval,
                                      M.ll = res$M$CI_M[1], M.ul = res$M$CI_M[2],
                                      A.ll = res$A$CI_A[1], A.ul = res$A$CI_A[2],
                                      phase.ll = res$phase$CI_phase[1], phase.ul = res$phase$CI_phase[2],
                                      Fstat = res$test$Fstat, sigma2 = res$test$sigma2, R2 = res$test$R2))
    }, parallel, cores)
    rhythmicity.tab = do.call(rbind.data.frame, rhythmicity.list)

  }else if(method == "permutation"){
    if((permutation.save!="NULL")&!dir.exists(permutation.save)){
      dir.create(file.path(permutation.save), recursive = TRUE)
      message(paste0("Directory has been created. Permutation results will be saved in ", permutation.save))
    }
    rhythmicity.tab = rhythmicity_permutation(X, t, period, nPermutation, permutation.save, permutation.file.label, parallel, cores)

  }else{
    stop("Please select method = 'OLS' or method = 'permutation'.")
  }

  rhythmicity.tab$qvalue = stats::p.adjust(rhythmicity.tab$pvalue, p.adjust.method)
  x$rhythm = rhythmicity.tab
  return(x)
}
