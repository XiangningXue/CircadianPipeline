#the least square cosinor package
#' Title Fit two-component cosinor with OLS for one gene
#'
#'This function fit the expression and time to this function \deqn{y = M1+A1(cos(2\pi/24*t+phase1))+X(M*+A*(cos(2\pi/24*t+phase*)))+\epsilon}, where \eqn{\epsilon \sim N(0, \sigma^2)}; X = 0 for group 1 and X = 1 for group 2.
#'For group 2 the formula can be converted to \deqn{y = M2+A2(cos(2\pi/24*t+phase2))+\epsilon}
#'
#' @param tod a numeric vector. Time of death (e.g. time of expression of the gene). length should be the same as y
#' @param y a numeric vector. Gene expression level.
#' @param group A vector indicating the group label. 0 for group 1 and 1 for group 2.
#' @param alpha a number between 0 to 1. The critical level for the confidence interval
#' @param period The length of the rhythmicity cycle. When it is 24 (default), the signal is circadian.
#'
#' @return {M1}{A list of OLS estimate of M1 and the lower and upper limit of M1 corresponding to the critical level alpha}
#' @return {A1}{A list of estimate of A1, sd of A1 and the lower and upper limit of A1 corresponding to the critical level alpha}
#' @return {phase1}{A list of estimate of phase1, sd of phase1 and the lower and upper limit of phase1 corresponding to the critical level alpha}
#' @return {M2}{A list of OLS estimate of M2 and the lower and upper limit of M2 corresponding to the critical level alpha}
#' @return {A2}{A list of estimate of A2, sd of A2 and the lower and upper limit of A2 corresponding to the critical level alpha}
#' @return {phase2}{A list of estimate of phase2, sd of phase2 and the lower and upper limit of phase2 corresponding to the critical level alpha}
#' @return {test}{A list of test result: F statistics, p-value and R2, \eqn{\sigma^2} (the variance of the error term)}
#'
#' @export
#' @examples
#' #Simulate a time series data with 20 time points
#'
#' x.time = runif(20, min = 0, max = 24)
#' m = 5; A = 3; phase = pi/4; sigma = 1
#' y = m+A*cos(2*pi/24*x.time+phase)+rnorm(20, 0, sigma)
#' est = one_cosinor_OLS(tod = x.time, y = y, alpha = 0.05, period = 24)
#'
two_cosinor_OLS = function(tod = time, y = y, group, alpha = 0.05, period = 24){

  #alpha is the critical level for equal tailed CI
  n = length(tod)
  x1 = cos(2*pi*tod/period)
  x2 = sin(2*pi*tod/period)
  x3 = group*x1
  x4 = group*x2

  mat.X = matrix(c(rep(1, n), x1, x2, group,x3, x4), ncol = 6, byrow = FALSE)
  mat.S = t(mat.X)%*%mat.X
  # mat.XX = t(mat.X)%*%mat.X#mat.XX = mat.S
  # mat.S = matrix(c(n, sum(x1), sum(x2),
  #                  sum(x1), sum(x1^2), sum(x1*x2),
  #                  sum(x2), sum(x1*x2), sum(x2^2)),
  #                nrow = 3, byrow = TRUE)
  #vec.d = c(sum(y), sum(y*x1), sum(y*x2))
  vec.d = t(mat.X)%*%y

  #estimates for group 1
  mat.S.inv = solve(mat.S)
  est = mat.S.inv%*%vec.d
  m.hat = est[1]
  beta1.hat = est[2]
  beta2.hat = est[3]

  # truth$phase[2]; truth$amplitude[2]; truth$M[502]
  m1.hat = m.hat
  A1.hat = sqrt(beta1.hat^2 + beta2.hat^2)
  phase1.res = get_phase(beta1.hat, beta2.hat)
  phase1.hat = phase1.res$phase
  peak1 <- (period-period*phase1.hat/(2*pi))%%period
  #if(peak > period/4*3) peak = peak - period #not used

  #estimates for group 2
  dm.hat = est[4]
  dbeta1.hat = est[5]
  dbeta2.hat = est[6]

  m2.hat = m1.hat + dm.hat
  A2.hat = sqrt((beta1.hat+dbeta1.hat)^2 + (beta2.hat+dbeta2.hat)^2)
  phase2.res = get_phase(beta1.hat+dbeta1.hat, beta2.hat+dbeta2.hat)
  phase2.hat = phase2.res$phase
  peak2 <- (period-period*phase2.hat/(2*pi))%%period

  #inference
  TSS = sum((y-mean(y))^2)
  yhat1 = m.hat + beta1.hat*x1+beta2.hat*x2 #prediction from partial model
  yhat2 = m.hat + beta1.hat*x1+beta2.hat*x2 + dm.hat*group + dbeta1.hat*x3+dbeta2.hat*x4
  RSS1 = sum((y-yhat1)^2)
  RSS2 = sum((y-yhat2)^2)
  MSS1 = TSS-RSS1
  MSS2 = TSS-RSS2
  Fstat = ((RSS1-RSS2)/3)/(RSS2/(n-6))
  pval = stats::pf(Fstat, 3, n-6, lower.tail = FALSE)
  sigma2.hat = RSS2/(n-6)
  sigma.hat = sqrt(sigma2.hat)

  CI.m1.hat.radius = calculate_CI.M(mat.S.inv, A.t = matrix(c(1, 0, 0, 0, 0, 0), nrow = 1),
                            r.full = 6, n = length(tod), alpha, sigma.hat = sigma.hat)
  CI.m2.hat.radius = calculate_CI.M(mat.S.inv, A.t = matrix(c(1, 0, 0, 1, 0, 0), nrow = 1),
                                    r.full = 6, n = length(tod), alpha, sigma.hat = sigma.hat)

  se.hat.group1 = calculate_CI_A.phase.Taylor(mat.S.inv, 1, phase1.hat, A1.hat, sigma.hat)
  se.hat.group2 = calculate_CI_A.phase.Taylor(mat.S.inv, 2, phase2.hat, A2.hat, sigma.hat)

  CI_A.phase.Scheffe.group1 = calculate_CI_A.phase.Scheffe(mat.S.inv, rbind(c(0, 1, 0, 0, 0, 0),
                                                                            c(0, 0, 1, 0, 0, 0)),
                                                           sigma2.hat,
                                                           est, r.full = 6, n, alpha)
  CI_A.phase.Scheffe.group2 = calculate_CI_A.phase.Scheffe(mat.S.inv, rbind(c(0, 1, 0, 0, 1, 0),
                                                                            c(0, 0, 1, 0, 0, 1)),
                                                           sigma2.hat,
                                                           est, r.full = 6, n, alpha)
  M.overlap = CI_overlap(c(m1.hat-CI.m1.hat.radius, m1.hat+CI.m1.hat.radius), c(m2.hat-CI.m2.hat.radius, m2.hat+CI.m2.hat.radius))
  A.overlap = CI_overlap(CI_A.phase.Scheffe.group1$CI_A, CI_A.phase.Scheffe.group2$CI_A)
  phase.overlap = CI_overlap(CI_A.phase.Scheffe.group1$CI_phase, CI_A.phase.Scheffe.group2$CI_phase)
  M.overlap_CI1_est2 = CI1_est2_overlap(c(m1.hat-CI.m1.hat.radius, m1.hat+CI.m1.hat.radius), m2.hat)
  A.overlap_CI1_est2 = CI1_est2_overlap(CI_A.phase.Scheffe.group1$CI_A, A2.hat)
  phase.overlap_CI1_est2 = CI1_est2_overlap(CI_A.phase.Scheffe.group1$CI_phase, phase2.hat)

  #output
  out = list(g1 = list(M = list(est = m1.hat,
                                CI_M = c(m1.hat-CI.m1.hat.radius, m1.hat+CI.m1.hat.radius)),
                       A = list(est = A1.hat,
                                sd = se.hat.group1$se.A.hat,
                                CI_A_delta = c(A1.hat-stats::qt(1-alpha/2, n-6)*se.hat.group1$se.A.hat, A1.hat+stats::qt(1-alpha/2, n-6)*se.hat.group1$se.A.hat),
                                CI_A = CI_A.phase.Scheffe.group1$CI_A), #conservative CI
                       phase = list(est = phase1.hat,
                                    sd = se.hat.group1$se.phase.hat,
                                    CI_phase_delta = c(phase1.hat-stats::qt(1-alpha/2, n-6)*se.hat.group1$se.phase.hat, phase1.hat+stats::qt(1-alpha/2, n-6)*se.hat.group1$se.phase.hat),
                                    #tan = phase.res$tan,
                                    #CI_tan = c(phi.lower.limit$tan, phi.upper.limit$tan),
                                    CI_phase = CI_A.phase.Scheffe.group1$CI_phase), #conservative CI
                       peak = peak1),
             g2 = list(M = list(est = m2.hat,
                                CI_M = c(m2.hat-CI.m2.hat.radius, m2.hat+CI.m2.hat.radius)),
                       A = list(est = A2.hat,
                                sd = se.hat.group2$se.A.hat,
                                CI_A_delta = c(A2.hat-stats::qt(1-alpha/2, n-6)*se.hat.group2$se.A.hat, A2.hat+stats::qt(1-alpha/2, n-6)*se.hat.group2$se.A.hat),
                                CI_A = CI_A.phase.Scheffe.group1$CI_A), #conservative CI
                       phase = list(est = phase2.hat,
                                    sd = se.hat.group2$se.phase.hat,
                                    CI_phase_delta = c(phase2.hat-stats::qt(1-alpha/2, n-6)*se.hat.group2$se.phase.hat, phase2.hat+stats::qt(1-alpha/2, n-6)*se.hat.group2$se.phase.hat),
                                    #tan = phase.res$tan,
                                    #CI_tan = c(phi.lower.limit$tan, phi.upper.limit$tan),
                                    CI_phase = CI_A.phase.Scheffe.group2$CI_phase), #conservative CI
                       peak = peak2),
             test = list(global.pval = pval,
                         M.ind = c(M.overlap, M.overlap_CI1_est2),
                         A.ind = c(A.overlap, A.overlap_CI1_est2),
                         phase.ind = c(phase.overlap, phase.overlap_CI1_est2)))
  return(out)

}

# B.mat = rbind(c(sum((x1-mean(x1))^2), sum((x1-mean(x1))*(x2-mean(x2)))),
#               c(sum((x1-mean(x1))*(x2-mean(x2))), sum((x2-mean(x2))^2)))
#
# #CI (derive conservative CI for phi)
# B11 = sum((x1-mean(x1))^2)
# B12 = sum((x1-mean(x1))*(x2-mean(x2)))
# B22 = sum((x2-mean(x2))^2)
# #this equals cov(cbind(x1, x2))*(n-1) which is close to t(cbind(x1, x2))%*%cbind(x1, x2)
# #test 1
# A.t = rbind(c(0, 1, 0),
#             c(0, 0, 1))
# xx1 = cbind(1, x1, x2)
# A.t%*%t(xx1)%*%xx1%*%t(A.t) #similar as t(cbind(x1, x2))%*%cbind(x1, x2)
# #then
# xx.inv = solve(t(xx1)%*%xx1)
# solve(A.t%*%xx.inv%*%t(A.t)) #still the same

CI_overlap = function(CI1 = c(1, 1.4), CI2 = c(1.5, 2)){
  ind1 = (CI1[1]>CI2[1])&(CI1[1]<CI2[2])
  ind2 = (CI1[2]>CI2[1])&(CI1[2]<CI2[2])
  return(!(ind1|ind2))
}

CI1_est2_overlap = function(CI1 = c(1, 1.4), est2 = 2){
  ind = (est2>CI1[1])&(est2<CI1[2])
  return(!ind)
}


