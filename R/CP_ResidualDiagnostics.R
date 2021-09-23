#' Title
#'
#' @param x the output fro CP_Rhythmicity
#' @param genes.diag A vector of genes that wish to do residual diagnostics. If NULL, top 10 most rhythmic genes will be selcted.
#' @param plot If TRUE, plots will be output. Otherwise, only ggplots objects are returned.
#' @param filename If NULL, no file will be saved. If input a valid filename, the residual plots will be saved with the filename.
#' @param height height for the saved file. Only used if filename is not NULL.
#' @param width width for the saved file. Only used if filename is not NULL.
#'
#' @return A list with one list for each gene.
#' @export
#'
#' @examples
#'
#' x1.time = runif(20, min = 0, max = 24)
#' m1 = rnorm(100, 5); A1 = rnorm(100, 3); phase1 = runif(100, min = 0, max = 2*pi); sigma1 = 1
#' noise.mat1 = matrix(rnorm(100*20, 0, sigma1), ncol = 20, nrow = 100)
#' signal.mat1 = t(sapply(1:100, function(a){m1[a]+A1[a]*cos(2*pi/24*x1.time+phase1[a])}))
#' x1 = list(data = as.data.frame(noise.mat1 + signal.mat1),
#'           tod = x1.time,
#'           label = paste("gene", seq_len(100)))
#' x1.Rhythm = CP_Rhythmicity(x1, parallel = FALSE)
#' CP_ResidualDiagnostics(x1.Rhythm)
#'
CP_ResidualDiagnostics = function(x = x1.rhythm, genes.diag = NULL, plot = TRUE,
                                  filename = NULL, height = 8, width = 8){
  if(is.null(genes.diag)){
    genes.diag = x$label[order(x$rhythm$pvalue)[1:10]]
  }

  tod = x$tod
  diag.list = lapply(1:length(genes.diag), function(a){
    a.gene = genes.diag[a]
    a.data = as.numeric(x$data[match(a.gene, x$label), ])
    a.est = x$rhythm[match(a.gene, x$label), ]
    a.fitted = fun.cosior(tod, a.est$P, a.est$A, a.est$phase, a.est$M)
    a.resid = a.data-a.fitted
    #residuals against fitted values
    df = data.frame(FittedValues = a.fitted, Residuals = a.resid)
    p1 = ggplot2::ggplot(df, ggplot2::aes(x = FittedValues, y = Residuals))+
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed")+
      ggplot2::xlab("Fitted values")+
      ggplot2::geom_point()+
      ggplot2::geom_smooth()+
      ggplot2::ggtitle("Residuals vs Fitted",  subtitle = a.gene)+
      ggplot2::theme_bw()

    #QQ plot
    n = length(tod)
    x1 = cos(2*pi*tod/a.est$P)
    x2 = sin(2*pi*tod/a.est$P)
    mat.X = matrix(c(rep(1, n), x1, x2), ncol = 3, byrow = FALSE)
    hat.X = mat.X%*%solve(t(mat.X)%*%mat.X)%*%t(mat.X)
    a.resid.std = a.resid/sqrt(a.est$sigma2*(1-diag(hat.X)))

    p2 = ggplot2::ggplot(df, ggplot2::aes(sample = a.resid.std))+
      ggplot2::geom_qq()+ggplot2::geom_qq_line(linetype = "dashed")+
      ggplot2::ylab("Standarized residuals")+
      ggplot2::xlab("Theoretical Quantiles")+
      ggplot2::ggtitle("Normal Q-Q",  subtitle = a.gene)+
      ggplot2::theme_bw()

    #scale-location
      #hat matrix

    df3 = data.frame(FittedValues = a.fitted, Residuals.std = a.resid.std)
    p3 = ggplot2::ggplot(df3, ggplot2::aes(x = FittedValues, y = sqrt(abs(Residuals.std))))+
      ggplot2::xlab("Fitted values")+
      ggplot2::ylab("sqrt(abs(Standard residuals))")+
      ggplot2::geom_point()+
      ggplot2::geom_smooth()+
      ggplot2::ggtitle("Scale-Location",  subtitle = a.gene)+
      ggplot2::theme_bw()

    df4 = data.frame(Residuals = a.resid, den =dnorm(a.resid, mean(a.resid), sd(a.resid)))
    p4 = ggplot2::ggplot(df4, ggplot2::aes(x = Residuals))+
      ggplot2::geom_histogram(ggplot2::aes(y = ..density..) ,bins=10)+
      ggplot2::geom_line(ggplot2::aes(x = Residuals, y = den), color = "red")+
      ggplot2::ggtitle("Residual histogram",  subtitle = a.gene)+
      ggplot2::theme_bw()

    test = shapiro.test(a.resid)
    return(list(p1, p2, p3, p4, test))
  })

  names(diag.list) = genes.diag

  pvals0 = lapply(diag.list, `[[`, 5)
  pvals = sapply(pvals0, `[[`, "p.value")
  if(all(pvals>0.05)){
    a.message = "All genes pass the Shapiro-Wilk normality test"
    message(a.message)
  }else{
    a.message = paste("Genes that do not pass the Shapiro-Wilk normality test:", paste0(names(diag.list)[pvals<0.05], collapse = ", "))
    message(a.message)
  }

  if(plot){
    suppressMessages(  lapply(diag.list, function(a){
      oask <- devAskNewPage(TRUE&dev.interactive()) #if in rmarkdown then will no longer ask. Otherwise hit to see the next plot
      on.exit(devAskNewPage(oask))
      gridExtra::grid.arrange(a[[1]], a[[2]], a[[3]], a[[4]], ncol = 2)
    }))
  }

  if(!is.null(filename)){
    pdf(filename, width, height)
    suppressMessages(  lapply(diag.list, function(a){
      gridExtra::grid.arrange(a[[1]], a[[2]], a[[3]], a[[4]], ncol = 2)
    }))
    dev.off()
  }
  #residuals against fitted values
  return(diag.list)
}

fun.cosior = function(x, P, A, phase, M){
  M+A*cos(2*pi/P*x+phase)
}


