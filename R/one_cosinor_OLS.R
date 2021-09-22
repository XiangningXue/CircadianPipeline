#the least square cosinor package
#' Title Fit one-component cosinor with OLS for one gene
#'
#' This function fit the expression and time to this function \deqn{y = M+A(cos(2\pi/24*t+phase))+\epsilon}, where \eqn{\epsilon \sim N(0, \sigma^2)}.
#'
#' @param tod a numeric vector. Time of death (e.g. time of expression of the gene). length should be the same as y
#' @param y a numeric vector. Gene expression level.
#' @param alpha a number between 0 to 1. The critical level for the confidence interval
#' @param period The length of the rhythmicity cycle. When it is 24 (default), the signal is circadian.
#' @param CI If return CI
#' @param CI.type "conservative": a convervative CI will be returned; "delta": CI calculated with sd derived with delta method and the plug in estimate will be returned.
#'
#' @return {M}{A list of OLS estimate of M and the lower and upper limit of M corresponding to the critical level alpha}
#' @return {A}{A list of estimate of A, sd of A and the lower and upper limit of A corresponding to the critical level alpha}
#' @return {phase}{A list of estimate of phase, sd of phase and the lower and upper limit of phase corresponding to the critical level alpha}
#' @return {peak}{The conversion from phase to peak is  \eqn{peak = (period-period*phase/(2*\pi))}.}
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
one_cosinor_OLS = function(tod = time, y = y, alpha = 0.05, period = 24, CI = TRUE, CI.type = "conservative"){

  #alpha is the critical level for equal tailed CI
  n = length(tod)
  x1 = cos(2*pi*tod/period)
  x2 = sin(2*pi*tod/period)
  # mat.X = matrix(c(rep(1, n), x1, x2), ncol = 3, byrow = FALSE)
  # mat.XX = t(mat.X)%*%mat.X#mat.XX = mat.S
  mat.S = matrix(c(n, sum(x1), sum(x2),
                   sum(x1), sum(x1^2), sum(x1*x2),
                   sum(x2), sum(x1*x2), sum(x2^2)),
                 nrow = 3, byrow = TRUE)
  vec.d = c(sum(y), sum(y*x1), sum(y*x2))

  mat.S.inv = solve(mat.S)
  est = mat.S.inv%*%vec.d
  m.hat = est[1]
  beta1.hat = est[2]
  beta2.hat = est[3]
  # truth$phase[2]; truth$amplitude[2]; truth$M[502]
  A.hat = sqrt(beta1.hat^2 + beta2.hat^2)

  phase.res = get_phase(beta1.hat, beta2.hat)
  phase.hat = phase.res$phase

  peak <- (period-period*phase.hat/(2*pi))%%period
  #if(peak > period/4*3) peak = peak - period #not used

  #inference
  TSS = sum((y-mean(y))^2)
  yhat = m.hat + beta1.hat*x1+beta2.hat*x2
  RSS = sum((y-yhat)^2)
  MSS = TSS-RSS
  Fstat = (MSS/2)/(RSS/(n-3))
  pval = stats::pf(Fstat, 2, n-3, lower.tail = FALSE)
  sigma2.hat = RSS/(n-3)
  sigma.hat = sqrt(sigma2.hat)

  if(CI){
    #CI (M and se for A and phi)
    CI.m.hat.radius = calculate_CI.M(mat.S.inv, A.t = matrix(c(1, 0, 0), nrow = 1),
                                      r.full = 3, n, alpha, sigma.hat)
    se.hat.A.phase = calculate_CI_A.phase.Taylor(mat.S.inv, 1, phase.hat, A.hat, sigma.hat)

    CI_A.phase.Scheffe = calculate_CI_A.phase.Scheffe(mat.S.inv, rbind(c(0, 1, 0),
                                                                       c(0, 0, 1)),
                                                      sigma2.hat,
                                                      est, r.full=3, n, alpha, "PlugIn") #The PlutIn CI is conservative enough.
    if(CI.type=="conservative"){
      a.CI = CI_A.phase.Scheffe
    }else{
      a.CI = list(CI_A = c(A.hat-stats::qt(1-alpha/2, n-3)*se.hat.A.phase$se.A.hat, A.hat+stats::qt(1-alpha/2, n-3)*se.hat.A.phase$se.A.hat),
                  CI_phase = c(phase.hat-stats::qt(1-alpha/2, n-3)*se.hat.A.phase$se.phase.hat, phase.hat+stats::qt(1-alpha/2, n-3)*se.hat.A.phase$se.phase.hat))
    }
    #output
    out = list(M = list(est = m.hat,
                        CI_M = c(m.hat-CI.m.hat.radius, m.hat+CI.m.hat.radius)),
               A = list(est = A.hat,
                        sd = se.hat.A.phase$se.A.hat,
                        CI_A = a.CI$CI_A), #conservative CI
               phase = list(est = phase.hat,
                            sd = se.hat.A.phase$se.phase.hat,
                            #tan = phase.res$tan,
                            #CI_tan = c(phi.lower.limit$tan, phi.upper.limit$tan),
                            CI_phase = a.CI$CI_phase), #conservative CI
               peak = peak,
               test = list(Fstat = Fstat,
                           pval = pval,
                           R2 = MSS/TSS,
                           sigma2 = sigma2.hat))
  }else{
    #output
    out = list(M = list(est = m.hat),
               A = list(est = A.hat), #conservative CI
               phase = list(est = phase.hat), #conservative CI
               peak = peak,
               test = list(R2 = MSS/TSS,
                           sigma2 = sigma2.hat))
  }
 return(out)
}

calculate_CI.M = function(XX.inv = mat.S.inv, A.t = matrix(c(1, 0, 0, 1, 0, 0), nrow = 1),
                          r.full = 6, n = length(tod), alpha = 0.05, sigma.hat = sigma.hat){
  CI.m.hat.radius = stats::qt(1-alpha/2, n-r.full)*sigma.hat*sqrt(A.t%*%XX.inv%*%t(A.t))
  return(CI.m.hat.radius)
}

calculate_CI_A.phase.Taylor = function(XX.inv = mat.S.inv, ind.group = 1, phase.hat = phase1.hat, A.hat = A1.hat,
                                       sigma.hat = sigma.hat){
  if(ind.group==1){
    var.beta1 = XX.inv[2, 2]; var.beta2 = XX.inv[3, 3]; var.beta1.beta2 = XX.inv[2, 3]
  }else if(ind.group==2){
    A.t = rbind(c(0, 1, 0, 0, 1, 0),
                c(0, 0, 1, 0, 0, 1))
    var.new = A.t%*%XX.inv%*%t(A.t)
    var.beta1 = var.new[1, 1]; var.beta2 = var.new[2, 2]; var.beta1.beta2 = var.new[1, 2]
  }
  se.A.hat = sigma.hat*sqrt(var.beta1*cos(phase.hat)^2
                            -2*var.beta1.beta2*sin(phase.hat)*cos(phase.hat)
                            +var.beta2*sin(phase.hat)^2)
  se.phase.hat = sigma.hat*sqrt(var.beta1*sin(phase.hat)^2
                                +2*var.beta1.beta2*sin(phase.hat)*cos(phase.hat)
                                +var.beta2*cos(phase.hat)^2)/A.hat
  return(list(se.A.hat = se.A.hat,
              se.phase.hat = se.phase.hat))

}

calculate_CI_A.phase.Scheffe = function(XX.inv = mat.S.inv, A.t = rbind(c(0, 1, 0, 0, 1, 0),
                                                                        c(0, 0, 1, 0, 0, 1)),
                                        sigma2.hat = sigma2.hat,
                                        est = est,r.full = 6, n = length(tod), alpha, CItype = "conservative"){
  #CItype = "conservative" or "PlugIn"
  # XX.inv = mat.S.inv;
  # A.t = rbind(c(0, 1, 0),
  #             c(0, 0, 1))
  # r.full=3

  q = Matrix::rankMatrix(A.t)[[1]]
  est2 = A.t%*%est
  beta1.hat = est2[1]
  beta2.hat = est2[2]
  B.inv = solve(A.t%*%XX.inv%*%t(A.t))
  B11 = B.inv[1, 1]; B12 = B.inv[1, 2]; B22 = B.inv[2, 2]
  R = q*sigma2.hat*stats::qf(1-alpha, q, n-r.full)

  C1 = -(B11*beta1.hat+B12*beta2.hat)/(B22*beta2.hat+B12*beta1.hat)
  C2 = -(R-2*B12*beta1.hat*beta2.hat-B11*beta1.hat^2-B22*beta2.hat^2)/(B22*beta2.hat+B12*beta1.hat)
  D1 = B22*C1^2+B11+2*B12*C1
  D2 = 2*B22*C1*C2+2*B12*C2-(B12*beta1.hat+B22*beta2.hat)*C1-B11*beta1.hat-B12*beta2.hat
  D3 = B22*C2^2-(B12*beta1.hat+B22*beta2.hat)*C2

  #calculate CI of phi
  #check if 0 is in ellipse
  zero.in.ellipse = B11*beta1.hat^2+2*B12*beta1.hat*beta2.hat+B22*beta2.hat^2 < R
  if(CItype == "conservative"){
    if(!zero.in.ellipse){
      delta.poly = D2^2-4*D1*D3
      if(delta.poly<0){
        phi.lower.limit = list(tan = -99, phase = -99)
        phi.upper.limit = list(tan = -99, phase = -99)
      }else{
        phi.beta1.roots = c((-D2-sqrt(delta.poly))/(2*D1), (-D2+sqrt(delta.poly))/(2*D1))
        phi.beta2.roots = C1*phi.beta1.roots+C2

        # phi.limit1 = get_phase(phi.beta1.roots[1], phi.beta2.roots[1])
        # phi.limit2 = get_phase(phi.beta1.roots[2], phi.beta2.roots[2])

        #change to adjusted phi limits
        phi.limits = get_phaseForCI(b1.x = beta1.hat, b2.x = beta2.hat,
                                    b1.r1 = phi.beta1.roots[1], b2.r1 = phi.beta2.roots[1],
                                    b1.r2 = phi.beta1.roots[2], b2.r2 = phi.beta2.roots[2])
        phi.lower.limit = phi.limits$phi.lower.limit
        phi.upper.limit = phi.limits$phi.upper.limit
      }
    }else{
      phi.lower.limit = list(tan = 99, phase = 99)
      phi.upper.limit = list(tan = 99, phase = 99)
    }

    #calculate CI of A
    #calculate normal ellipse parameters
    ellipse.parameters = solve.ellipse.parameters(a = B11, b = B22, c = 2*B12,
                                                  d = -2*(B11*beta1.hat+B12*beta2.hat),
                                                  e = -2*(B12*beta1.hat+B22*beta2.hat),
                                                  f = B11*beta1.hat^2+B22*beta2.hat^2+2*B12*beta1.hat*beta2.hat-q*sigma2.hat*stats::qf(1-alpha, q, n-r.full))
    angle.point.newOrigin.major =
      get.angle_point.newOrigin.major(x0 = ellipse.parameters$x0,
                                      y0 = ellipse.parameters$y0,
                                      theta.rotate = ellipse.parameters$theta.rotate)
    r.point.to.center = sqrt(ellipse.parameters$x0^2+ellipse.parameters$y0^2)
    x.new = r.point.to.center*cos(angle.point.newOrigin.major)
    y.new = r.point.to.center*sin(angle.point.newOrigin.major)

    #check the position of the new origin to the ellipse
    if(x.new==0&y.new==0){
      A.limit1 = ellipse.parameters$minor
      A.limit2 = ellipse.parameters$major
    }else if(x.new==0){
      A.limit1 = abs(ellipse.parameters$minor-y.new)
      A.limit2 = ellipse.parameters$minor+y.new
    }else if(y.new==0){
      A.limit1 = abs(ellipse.parameters$major-x.new)
      A.limit2 = ellipse.parameters$major+x.new
    }else if(x.new!=0&y.new!=0){
      x.new = abs(x.new)
      y.new = abs(y.new)
      fun.t = function(t){
        (ellipse.parameters$major*x.new/(t+ellipse.parameters$major^2))^2+
          (ellipse.parameters$minor*y.new/(t+ellipse.parameters$minor^2))^2-1
      }

      # root.upper = 50
      # while(fun.t(-ellipse.parameters$minor^2+0.0001)*fun.t(root.upper)>0){
      #   root.upper = root.upper+50
      # }
      #    root1 = uniroot(fun.t, c(-ellipse.parameters$minor^2, root.upper),extendInt="downX")$root
      root1 = stats::uniroot(fun.t, c(-ellipse.parameters$minor^2, -ellipse.parameters$minor^2+50),extendInt="downX")$root
      x.root1 = ellipse.parameters$major^2*x.new/(root1+ellipse.parameters$major^2)
      y.root1 = ellipse.parameters$minor^2*y.new/(root1+ellipse.parameters$minor^2)
      dmin = sqrt((x.new-x.root1)^2+(y.new-y.root1)^2)

      # root.lower = -50
      # while(fun.t(-ellipse.parameters$major^2-0.0001)*fun.t(root.lower)>0){
      #   root.lower = root.lower-50
      # }
      #root2 = uniroot(fun.t, c(root.lower, -ellipse.parameters$major^2-0.0001), extendInt="yes")$root
      root2 = stats::uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="upX")$root
      # root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="yes")$root
      # root2 = uniroot(fun.t, c(-ellipse.parameters$major^2-50, -ellipse.parameters$major^2), extendInt="no")$root
      x.root2 = ellipse.parameters$major^2*x.new/(root2+ellipse.parameters$major^2)
      y.root2 = ellipse.parameters$minor^2*y.new/(root2+ellipse.parameters$minor^2)
      dmax = sqrt((x.new-x.root2)^2+(y.new-y.root2)^2)

      A.limit1 = min(dmin, dmax)
      A.limit2 = max(dmin, dmax)
    }
  }else if(CItype=="PlugIn"){
    if(!zero.in.ellipse){
      delta.poly = D2^2-4*D1*D3
      if(delta.poly<0){
        phi.lower.limit = list(tan = -99, phase = -99)
        phi.upper.limit = list(tan = -99, phase = -99)
      }else{
        #get the roots through the rootSolve package
        #install.packages("rootSolve")
        phi.beta1.roots = c((-D2-sqrt(delta.poly))/(2*D1), (-D2+sqrt(delta.poly))/(2*D1))
        phi.beta2.roots = C1*phi.beta1.roots+C2

        # model = function(x){
        #   beta1 = x[1]; beta2 = x[2]
        #   F1 = B11*(beta1-beta1.hat)^2+2*B12*(beta1-beta1.hat)*(beta2-beta2.hat)+B22*(beta2-beta2.hat)^2-R
        #   F2 = beta1^2+beta2^2-beta1.hat^2-beta2.hat^2
        #   return(c(F1 = F1, F2 = F2))
        # }
        #
        # start.points = list(c(phi.beta1.roots[1], phi.beta2.roots[1]),
        #                     c(phi.beta1.roots[2], phi.beta2.roots[2]),
        #                     c(beta1.hat/2, beta2.hat/2),
        #                     c(beta1.hat/2, beta2.hat*2),
        #                     c(beta1.hat*2, beta2.hat/2),
        #                     c(beta1.hat*2, beta2.hat*2))
        #
        # solutions.A.PlugIn = lapply(start.points, function(a){
        #   res = tryCatch(rootSolve::multiroot(f=model, start = a , maxiter = 100)$root, warning=function(cnd){NA})
        #   return(res)
        # })
        #
        # dist.solutions = as.matrix(dist(do.call(rbind, solutions.A.PlugIn)))
        # solution1.A.PlugIn = solutions.A.PlugIn[[1]]
        # solution2.A.PlugIn = solutions.A.PlugIn[[which(round(dist.solutions[, 1], 1)!=0)[1]]]
        A.hat2 = sqrt(beta1.hat^2+beta2.hat^2)

        fun = function(phi){
          B11*A.hat2^2*cos(phi)^2+B22*A.hat2^2*sin(phi)^2+2*B12*A.hat2^2*sin(phi)*cos(phi)-
            (2*B11*beta1.hat*A.hat2+2*B12*beta2.hat*A.hat2)*cos(phi)-
            (2*B12*beta1.hat*A.hat2+2*B22*beta2.hat*A.hat2)*sin(phi)+
            B11*beta1.hat^2+2*B12*beta1.hat*beta2.hat+B22*beta2.hat^2-R
        }
        all.solutions = rootSolve::uniroot.all(fun, c(0, 2*pi))

        solution1.A.PlugIn = c(A.hat2*cos(all.solutions[1]), A.hat2*sin(all.solutions[1]))
        solution2.A.PlugIn = c(A.hat2*cos(all.solutions[2]), A.hat2*sin(all.solutions[2]))


        #change to adjusted phi limits
        phi.limits = get_phaseForCI(b1.x = beta1.hat, b2.x = beta2.hat,
                                    b1.r1 = solution1.A.PlugIn[1], b2.r1 = solution1.A.PlugIn[2],
                                    b1.r2 = solution2.A.PlugIn[1], b2.r2 = solution2.A.PlugIn[2])
        phi.lower.limit = phi.limits$phi.lower.limit
        phi.upper.limit = phi.limits$phi.upper.limit
      }
    }else{
      phi.lower.limit = list(tan = 99, phase = 99)
      phi.upper.limit = list(tan = 99, phase = 99)
    }

    #calculate CI of A
    #calculate normal ellipse parameters
    ellipse.parameters = solve.ellipse.parameters(a = B11, b = B22, c = 2*B12,
                                                  d = -2*(B11*beta1.hat+B12*beta2.hat),
                                                  e = -2*(B12*beta1.hat+B22*beta2.hat),
                                                  f = B11*beta1.hat^2+B22*beta2.hat^2+2*B12*beta1.hat*beta2.hat-q*sigma2.hat*stats::qf(1-alpha, q, n-r.full))
    angle.point.newOrigin.major =
      get.angle_point.newOrigin.major(x0 = ellipse.parameters$x0,
                                      y0 = ellipse.parameters$y0,
                                      theta.rotate = ellipse.parameters$theta.rotate)
    r.point.to.center = sqrt(ellipse.parameters$x0^2+ellipse.parameters$y0^2)
    x.new = r.point.to.center*cos(angle.point.newOrigin.major)
    y.new = r.point.to.center*sin(angle.point.newOrigin.major)

    #check the position of the new origin to the ellipse
    if(x.new==0&y.new==0){
      A.limit1 = ellipse.parameters$minor
      A.limit2 = ellipse.parameters$major
    }else if(x.new==0){
      A.limit1 = abs(ellipse.parameters$minor-y.new)
      A.limit2 = ellipse.parameters$minor+y.new
    }else if(y.new==0){
      A.limit1 = abs(ellipse.parameters$major-x.new)
      A.limit2 = ellipse.parameters$major+x.new
    }else if(x.new!=0&y.new!=0){
      x.new = abs(x.new)
      y.new = abs(y.new)
      r1.xx2 = ellipse.parameters$major^2*ellipse.parameters$minor^2/(ellipse.parameters$major^2+ellipse.parameters$minor^2*ellipse.parameters$tan.rotate^2)
      r1 = sqrt(r1.xx2+r1.xx2*ellipse.parameters$tan.rotate^2)
      dmin = sqrt(x.new^2+y.new^2)-r1
      dmax = sqrt(x.new^2+y.new^2)+r1

      A.limit1 = min(dmin, dmax)
      A.limit2 = max(dmin, dmax)
    }
  }

  if(zero.in.ellipse){
    A.limit1 = 0
  }

  return(list(CI_phase = c(phi.lower.limit$phase, phi.upper.limit$phase),
              CI_A = c(A.limit1, A.limit2)))

}

get_phase = function(b1.x = beta1.hat, b2.x = beta2.hat){
  ph.x = atan(-b2.x/b1.x)
  #adjust ph.x
  if(b2.x>0){
    if(ph.x<0){
      ph.x = ph.x+2*pi
    }else if (ph.x>0){
      ph.x = ph.x+pi
    }
  }else if(b2.x<0){
    if(ph.x<0){
      ph.x = ph.x+pi
    }
  }else{
    ph.x = 88#88 means one the the beta estimate is 0. I did not account for such senerio because it is rare and complicated
  }
  return(list(phase = ph.x,
              tan = -b2.x/b1.x))
}


get_phaseForCI = function(b1.x = beta1.hat, b2.x = beta2.hat,
                          b1.r1 = phi.beta1.roots[1], b2.r1 = phi.beta2.roots[1],
                          b1.r2 = phi.beta1.roots[2], b2.r2 = phi.beta2.roots[2]){
  # b1.x = beta1.hat; b2.x = beta2.hat;
  # b1.r1 = phi.beta1.roots[1]; b2.r1 = phi.beta2.roots[1];
  # b1.r2 = phi.beta1.roots[2]; b2.r2 = phi.beta2.roots[2]

  #b1.r1 is the beta1 estimate of root 1.
  #step 1: find which root is in the same quadrant as OLS estimate b1.x, b2.x

  x.quad = get_quad(b1.q = b1.x, b2.q = b2.x)
  phase.res = get_phase(b1.x, b2.x)
  r1.quad = get_quad(b1.q = b1.r1, b2.q = b2.r1)
  r2.quad = get_quad(b1.q = b1.r2, b2.q = b2.r2)

  ##the easy scenario: three are in the same quadrant
  if(r1.quad==x.quad&r2.quad==x.quad){
    phi.limit1 = get_phase(b1.r1, b2.r1)
    phi.limit2 = get_phase(b1.r2, b2.r2)
    if(phi.limit1$phase<phi.limit2$phase){
      phi.lower.limit = phi.limit1
      phi.upper.limit = phi.limit2
    }else{
      phi.lower.limit = phi.limit2
      phi.upper.limit = phi.limit1
    }
  }else if(r1.quad == x.quad|r2.quad==x.quad){
    #more complicated scenario: one of the root is in the other quadrant
    if(r1.quad == x.quad){
      b1.s1 = b1.r1; b2.s1 = b2.r1 #s1 stands for the root that is in the same quadrant
      b1.s2 = b1.r2; b2.s2 = b2.r2
      quad.s1 = r1.quad; quad.s2 = r2.quad
    }else{
      b1.s1 = b1.r2; b2.s1 = b2.r2 #s1 stands for the root that is in the same quadrant
      b1.s2 = b1.r1; b2.s2 = b2.r1
      quad.s1 = r2.quad; quad.s2 = r1.quad
    }
    phi.s1 = get_phase(b1.s1, b2.s1)
    if(phi.s1$phase<phase.res$phase){
      phi.lower.limit = phi.s1
      s2.type = "s1 lower, s2 upper"
      phi.upper.limit = get_phaseForCI_QuadTab(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type,
                                               b1.s2.q = b1.s2, b2.s2.q = b2.s2)
    }else{
      phi.upper.limit = phi.s1
      s2.type = "s1 upper, s1 lower"
      phi.lower.limit = get_phaseForCI_QuadTab(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type,
                                               b1.s2.q = b1.s2, b2.s2.q = b2.s2)
    }
  }else{
    #this is the scenario that the ellipse covers three quadrants
    phi.limits = get_phaseForCI_QuadTab2(x.quad, r1.quad, r2.quad,
                                         b1.t1 = b1.r1, b2.t1 = b2.r1,
                                         b1.t2 = b1.r2, b2.t2 = b2.r2)
    phi.lower.limit = phi.limits$phi.lower.limit
    phi.upper.limit = phi.limits$phi.upper.limit
  }

  return(list(phi.lower.limit = phi.lower.limit,
              phi.upper.limit = phi.upper.limit))

}

#make a table of the quadrants
get_quad = function(b1.q = b1.x, b2.q = b2.x){
  quad.table = as.data.frame(
    list(beta1.gt.0 = c(TRUE, TRUE, FALSE, FALSE),
         beta2.gt.0 = c(TRUE, FALSE, TRUE, FALSE),
         quad = c(1, 4, 2, 3))
  )
  quad = quad.table[(quad.table$beta1.gt.0==(b1.q>0))&(quad.table$beta2.gt.0==(b2.q>0)),
                    "quad"]
  return(quad)
}

# s1.quad = quad.s1; s2.quad = quad.s2; s1s2 = s2.type;
# b1.s2.q = b1.s2; b2.s2.q = b2.s2
get_phaseForCI_QuadTab = function(s1.quad = quad.s1, s2.quad = quad.s2, s1s2 = s2.type,
                                  b1.s2.q = b1.s2, b2.s2.q = b2.s2){
  # s1.quad = quad.s1; s2.quad = quad.s2; s1s2 = s2.type;
  # b1.s2.q = b1.s2; b2.s2.q = b2.s2
  quad.table2 = as.data.frame(
    list(s1.quad = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4),
         s2.quad = c(2, 2, 3, 3, 4, 4, 1, 1, 3, 3, 4, 4, 1, 1, 2, 2, 4, 4, 1, 1, 2, 2, 3, 3),
         s1s2 = rep(c("s1 lower, s2 upper", "s1 upper, s1 lower"), 12),
         s2.low = c(3, 1, 5/2, 1/2, 2, 0, 3/2, -1/2, 5/2, 1/2, 2, 0, 3/2, -1/2, 1, -1, 2, 0, 3/2, -1/2, 1, -1, 1/2, -3/2),
         s2.up  = c(7/2, 3/2, 3, 1, 5/2, 1/2, 2, 0, 3, 1, 5/2, 1/2, 2, 0, 3/2, -1/2, 5/2, 1/2, 2, 0, 3/2, -1/2, 1, -1))
  )

  s2.low = quad.table2[quad.table2$s1.quad==s1.quad&quad.table2$s2.quad==s2.quad&quad.table2$s1s2==s1s2,
                       "s2.low"]
  s2.up  = quad.table2[quad.table2$s1.quad==s1.quad&quad.table2$s2.quad==s2.quad&quad.table2$s1s2==s1s2,
                       "s2.up"]

  phi.s2 = atan(-b2.s2.q/b1.s2.q)
  if(phi.s2<s2.low*pi){
    while(phi.s2<s2.low*pi){
      phi.s2 = phi.s2+pi
    }
  }else if(phi.s2>s2.up*pi){
    while(phi.s2>s2.up*pi){
      phi.s2 = phi.s2-pi
    }
  }

  # if((phi.s2<(s2.up*pi))&(phi.s2>(s2.low*pi))){
  #   #print("#S2 is in interval")
  # }else{
  #   #print("S2 is not in the int!!!!!!!!!!!!!!!!")
  # }
  return(list(phase = phi.s2,
              tan = -b2.s2.q/b1.s2.q))
}

get_phaseForCI_QuadTab2 = function(x.quad, r1.quad, r2.quad,
                                   b1.t1 = b1.r1, b2.t1 = b2.r1,
                                   b1.t2 = b1.r2, b2.t2 = b2.r2){
  #t1 stands for three-quadrants, solution 1
  #This function is for when the ellipse covers three quadrants
  quad.table3 = as.data.frame(
    list(x.quad = c(1, 1, 2, 2, 3, 3, 4, 4),
         t.quad = c(2, 4, 3, 1, 4, 2, 1, 3),
         t.low = c(1, 2, 1/2, 3/2, 0, 1, -1/2, 1/2),
         t.up = c(3/2, 5/2, 1, 2, 1/2, 3/2, 0, 1))
  )
  t1.low = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r1.quad, "t.low"]
  t1.up = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r1.quad, "t.up"]
  t2.low = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r2.quad, "t.low"]
  t2.up = quad.table3[quad.table3$x.quad==x.quad&quad.table3$t.quad==r2.quad, "t.up"]

  phi.t1 = atan(-b2.t1/b1.t1)
  if(phi.t1<t1.low*pi){
    while(phi.t1<t1.low*pi){
      phi.t1 = phi.t1+pi
    }
  }else if(phi.t1>t1.up*pi){
    while(phi.t1>t1.up*pi){
      phi.t1 = phi.t1-pi
    }
  }

  # if((phi.t1<(t1.up*pi))&(phi.t1>(t1.low*pi))){
  #  # print("#t1 is in interval")
  # }else{
  #  # print("t1 is not in the int!!!!!!!!!!!!!!!!")
  # }

  phi.t2 = atan(-b2.t2/b1.t2)
  if(phi.t2<t2.low*pi){
    while(phi.t2<t2.low*pi){
      phi.t2 = phi.t2+pi
    }
  }else if(phi.t2>t2.up*pi){
    while(phi.t2>t2.up*pi){
      phi.t2 = phi.t2-pi
    }
  }

  # if((phi.t2<(t2.up*pi))&(phi.t2>(t2.low*pi))){
  #   # print("#t2 is in interval")
  # }else{
  #   # print("t2 is not in the int!!!!!!!!!!!!!!!!")
  # }

  if(phi.t1<phi.t2){
    phi.lower.limit = list(phase = phi.t1,
                           tan = -b2.t1/b1.t1)
    phi.upper.limit = list(phase = phi.t2,
                           tan = -b2.t2/b1.t2)
  }else{
    phi.lower.limit = list(phase = phi.t2,
                           tan = -b2.t2/b1.t2)
    phi.upper.limit = list(phase = phi.t1,
                           tan = -b2.t1/b1.t1)
  }
  return(list(phi.lower.limit = phi.lower.limit,
              phi.upper.limit = phi.upper.limit))
}

solve.ellipse.parameters = function(a = a, b = b, c = c, d = d, e = e, f = f){
  #a * x ^ 2 + b * y ^ 2 + c * x * y + d * x + e * y + f = 0
  delta1 = c^2 -4*a*b
  if(delta1>=0){
    return("Not a proper epplise")
  }else{
    #the transformation function is from wikepedia
    denominator.factor1 = 2*(a*e^2+b*d^2-c*d*e+delta1*f)
    major.minor.square.diff = sqrt(((a-b)^2+c^2))
    denominator.factor2 = c(a+b+major.minor.square.diff, a+b-major.minor.square.diff)
    major.minor = -sqrt(denominator.factor1*denominator.factor2)/delta1
    x0 = (2*b*d-c*e)/delta1
    y0 = (2*a*e-c*d)/delta1
    if(c==0){
      if(a<b){
        theta = 0
        tan.rotate = 0
      }else{
        theta = pi/2
        tan.rotate = Inf
      }
    }else{
      tan.rotate = (b-a-major.minor.square.diff)/c
      theta = atan(tan.rotate)
    }
    return(list(major = major.minor[1],
                minor = major.minor[2],
                x0 = x0,
                y0 = y0,
                theta.rotate = theta,
                tan.rotate = tan.rotate))
  }
}

get.angle_point.newOrigin.major = function(x0 = ellipse.parameters$x0,
                                           y0 = ellipse.parameters$y0,
                                           theta.rotate = ellipse.parameters$theta.rotate){
  #if both x0=0 and y0=0, then the distance is already there:
  #min distance = b; max distance = a
  if(x0==0&y0==0){
    return("ellipse is centered at (0, 0)")
  }else if(y0==0){
    angle = theta.rotate
    return(angle)
  }else if(x0==0){
    angle = pi/2-theta.rotate
    return(angle)
  }else{
    tan.center.0.x_axis = y0/x0
    theta.center.0.x_axis = atan(tan.center.0.x_axis)
    if(theta.center.0.x_axis*theta.rotate<0){
      angle = abs(theta.center.0.x_axis)+abs(theta.rotate)
      return(angle)
    }else if(theta.center.0.x_axis*theta.rotate>0){
      angle = abs(theta.center.0.x_axis-theta.rotate)
      return(angle)
    }
  }
}




