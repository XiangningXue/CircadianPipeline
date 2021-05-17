##' Fit cosinor function
##'
##' This function fit the expression and time to this function \deqn{y = M+A(cos(2\pi/24*t+phase))+\epsilon}, where no distribution assumption is made on \eqn{\epsilon}.
##' @title Fit cosinor function using optimization
##' @param tod a numeric vector. Time of death (e.g. time of expression of the gene). length should be the same as y
##' @param y a numeric vector. Gene expression level.
##' @param period period The length of the rhythmicity cycle. When it is 24 (default), the signal is circadian.
##' @param parStart Initial value for optimization purpose.
##' @return A list of A, phase, M, peak, B1, B2, SST, SSE, R2.
##' Formula 1: \eqn{y = M + A * cos(2\pi/period * t+phase)}
##' Formula 2: \eqn{y = M + B1 * cos(2\pi/period * t) + B2 * cos(2*\pi/period * t) }
##' \item{A}{Amplitude based on formula 1.}
##' \item{phase}{Phase based on formula 1, phase is restricted within (0, 2\eqn{\pi}).}
##' \item{M}{Basal level (vertical shift) based on formula 1 or on formula 2.}
##' \item{B1}{B1 based on formula 2.}
##' \item{B2}{B2 based on formula 2.}
##' \item{tss}{Total sum of square.}
##' \item{rss}{Residual sum of square, SSE/n is the MLE of the variance sigma square.}
##' \item{R2}{Pseudo R2 defined as (tss - rss)/tss.}
##'
##' @import minpack.lm
##' @export
##' @examples
##' set.seed(32608)
##' n <- 50
##' tt <- runif(n,0,24)
##' Amp <- 2
##' Phase <- 6
##' Offset <- 3
##' yy <- Amp * cos(2*pi/24 * tt + Phase) + Offset + rnorm(n,0,1)
##' fitCosCurve(tt, yy)

fitCosCurve <- function(tt, yy, period = 24, parStart = list(amp=3,phase=0, offset=0)){

  getPred <- function(parS, tt) {
    parS$amp * cos(2*pi/period * tt + parS$phase) + parS$offset
  }

  residFun <- function(p, yy, tt) yy - getPred(p,tt)

  nls.out <- nls.lm(par=parStart, fn = residFun, yy = yy,	tt = tt)

  apar <- nls.out$par

  amp0 <- apar$amp
  asign <- sign(amp0)
  ## restrict amp > 0
  amp <- amp0 * asign

  phase0 <- apar$phase
  #phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% period

  phase <- (phase0 + ifelse(asign==1,0,pi)) %% (2*pi)
  offset <- apar$offset

  peak <- (period-period*phase/(2*pi))%%period
  if(peak > period/4*3) peak = peak - period

  B1 <- amp0 * cos(2*pi/period)*cos(phase0)
  B2 <- -1*amp0 * sin(2*pi/period)*sin(phase0)

  rss <- sum(nls.out$fvec^2)
  tss <- sum((yy - mean(yy))^2)
  R2 <- 1 - rss/tss

  if(F){
    amp <- apar$amp
    phase <- apar$phase
    offset <- apar$offset
  }

  res <- list(A=amp, phase=phase, M=offset, peak=peak, R2=R2, sigma2 =  rss/(n-3), RSS=rss)
  res
}


