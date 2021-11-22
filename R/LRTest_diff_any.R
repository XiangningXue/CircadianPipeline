
# set.seed(32608)
# n <- 50
# tt1 <- runif(n,0,24)
# Amp1 <- 3
# Phase1 <- 5
# Offset1 <- 3
# yy1 <- Amp1 * sin(2*pi/24 * (tt1 + Phase1)) + Offset1 + rnorm(n,0,1)
# tt2 <- runif(n,0,24)
# Amp2 <- 3
# Phase2 <- 5
# Offset2 <- 3
# yy2 <- Amp2 * sin(2*pi/24 * (tt2 + Phase2)) + Offset2 + rnorm(n,0,1)
# LRTest_diff_any(tt1, yy1, tt2, yy2)
# period = 24

LRTest_diff_any = function(tt1, yy1, tt2, yy2, period = 24, FN = TRUE){


  n1 <- length(tt1)
  stopifnot(n1 == length(yy1))
  n2 <- length(tt2)
  stopifnot(length(tt2) == length(yy2))

  w <- 2*pi/period

  fit1 = one_cosinor_OLS(tt1, yy1, period, CI = FALSE)
  fit2 = one_cosinor_OLS(tt2, yy2, period, CI = FALSE)

  A1 <- fit1$A$est
  A2 <- fit2$A$est

  phase1 <- fit1$phase$est
  phase2 <- fit2$phase$est

  basal1 <- fit1$M$est
  basal2 <- fit2$M$est

  sigma2_1 <- fit1$test$R2
  sigma2_2 <- fit2$test$R2

  theta1 <- 1/sigma2_1
  theta2 <- 1/sigma2_2
#
#   p1 <- c(A1, phase1, basal1, theta1)
#   p2 <- c(A2, phase2, basal2, theta2)
#
#   x_Ha <- c(p1, p2)

  yhat1 = basal1+A1*cos(w*tt1+phase1)
  ll1_a <- log(theta1)/2
  ll1_b <- (yy1 - yhat1)^2 * theta1 / 2
  ll1 <- ll1_a - ll1_b

  yhat2 = basal2+A2*cos(w*tt2+phase2)
  ll2_a <- log(theta2)/2
  ll2_b <- (yy2 - yhat2)^2 * theta2 / 2
  ll2 <- ll2_a - ll2_b

  ll_Ha = sum(ll1) + sum(ll2)

  tt_H0 = c(tt1, tt2); yy_H0 = c(yy1, yy2)

  fit_H0 = one_cosinor_OLS(tt_H0, yy_H0, period, CI = FALSE)

  A_H0 <- fit_H0$A$est

  phase_H0 <- fit_H0$phase$est

  basal_H0 <- fit_H0$M$est

  sigma2_H0 <- fit_H0$test$R2
  theta_H0 <- 1/sigma2_H0

  yhat_H0 = basal_H0+A_H0*cos(w*tt_H0+phase_H0)
  ll_H0_a <- log(theta_H0)/2
  ll_H0_b <- (yy_H0 - yhat_H0)^2 * theta_H0 / 2
  ll_H0_0 <- ll_H0_a - ll_H0_b

  ll_H0 = sum(ll_H0_0)

  LR_stat <- -2*(ll_H0-ll_Ha)

  dfdiff <- 4
  if(!FN){
    pvalue <- pchisq(LR_stat,dfdiff,lower.tail = F)
  } else if(FN){
    r <- 4
    k <- 8
    n <- n1+n2
    Fstat <- (exp(LR_stat/n) - 1) * (n-k) / r
    pvalue <- pf(Fstat,df1 = r, df2 = n-k, lower.tail = F)
  } else{
    stop("FN has to be TRUE or FALSE")
  }

  res <- list(A1=A1, A2 = A2, A_c = A_H0,
              phase_1=phase1, phase_2=phase2, phase_c=phase_H0,
              M1 = basal1, M2 = basal2, M_c = basal_H0,
              sigma2_1 = sigma2_1, sigma2_2 = sigma2_2, sigma2_H0 = sigma2_H0,
              l0=ll_H0,
              la=ll_Ha,
              #df = dfdiff,
              stat=LR_stat,
              pvalue=pvalue)
  return(res)
}
