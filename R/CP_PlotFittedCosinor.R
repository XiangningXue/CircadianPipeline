#' Title
#'
#' @param x1 output of CP_Rhythmicity(). A list of data, tod, label and rhythmicity parameter estimtggplot2::aes
#' @param x2 output of CP_Rhythmicity(). If NULL, then only one plot will be made
#' @param genes.plot labels of genes that wish to be plotted. Labels must be consistent with those given in x1 or x1 and x2. If NULL, the top 10 most rhythmic genes from x1 will be used
#' @param category1 category of x1 if there are subgroups of samples. Must be a vector whose length equals number of samples in x1.
#' @param category2 category of x2 if there are subgroups of samples. Must be a vector whose length equals number of samples in x2.
#' @param category1.label will be used as legend title
#' @param category2.label will be used as legend title
#' @param Info1 used in title of the plot
#' @param Info2 used in title of the plot
#' @param height height of the output plot
#' @param width width of the output plot
#' @param filename filename for image output. If NULL, only ggplot object will be returned.
#'
#' @return A list which contains: a list of ggplot2 objects if x2 is null; two lists if x2 is not null.
#' @export
#'
#' @examples
#' #Simulate two time series data with 20 time points and 100 genes
#'
#' x1.time = runif(20, min = 0, max = 24)
#' m1 = rnorm(100, 5); A1 = rnorm(100, 3); phase1 = runif(100, min = 0, max = 2*pi); sigma1 = 1
#' noise.mat1 = matrix(rnorm(100*20, 0, sigma1), ncol = 20, nrow = 100)
#' signal.mat1 = t(sapply(1:100, function(a){m1[a]+A1[a]*cos(2*pi/24*x1.time+phase1[a])}))
#' x1 = list(data = as.data.frame(noise.mat1 + signal.mat1),
#' tod = x1.time,
#' label = paste("gene", seq_len(100)))
#' x1.rhythm = CP_Rhythmicity(x1, parallel = FALSE)
#'
#' x2.time = runif(20, min = 0, max = 24)
#' m2 = m1 + 2; A2 = A1+1; phase2 = phase1+pi/4; sigma = 1
#' noise.mat2 = matrix(rnorm(100*20, 0, sigma1), ncol = 20, nrow = 100)
#' signal.mat2 = t(sapply(1:100, function(a){m2[a]+A2[a]*cos(2*pi/24*x2.time+phase2[a])}))
#' x2 = list(data = as.data.frame(noise.mat2 + signal.mat2),
#' tod = x2.time,
#' label = paste("gene", seq_len(100)))
#' x2.rhythm = CP_Rhythmicity(x2, parallel = FALSE)
#'
#' CP_PlotFittedCosinor(x1.rhythm, x2.rhythm, genes.plot = NULL, Info1 = "data1", Info2 = "data2")
CP_PlotFittedCosinor = function(x1 = x1.rhythm, x2 = NULL, genes.plot = NULL,
                                Info1 = "data1", Info2 = "data2",
                                category1 = NULL, category1.label = NULL, category2 = NULL, category2.label = NULL,
                                filename = NULL, height = 8, width = 8){

  if((!is.data.frame(x1$data))|(!is.data.frame(x1$rhythm))){
    stop("x1: data and rhythmicity parameters must be dataframes")
  }
  if(length(x1$data)!=length(x1$tod)){
    stop("x1: Number of samples in data does not match that in tod.")
  }
  if(nrow(x1$data)!=length(x1$label)){
    stop("x1: Number of labels does not match number of feature in data.")
  }
  if(!is.null(category1)){
    if(ncol(x1$data)!=length(category1)){
      stop("length of category 1 does not equal number of samples in x1")
    }
  }
  if(is.null(category1.label)){
    category1.label = "category1"
  }
  if(is.null(genes.plot)){
    genes.plot = x1$label[order(x1$rhythm$pvalue)[1:10]]
  }
  p1.list = lapply(1:length(genes.plot), function(a){
    gene1 = genes.plot[a]
    tod1 = x1$tod
    expr1 = as.numeric(x1$data[match(gene1, x1$label),])
    apar1 = x1$rhythm[match(gene1, x1$label),]
    myylim = c(min(expr1), max(expr1))
    p1 = circadianDrawing_one(tod1, expr1, apar1, gene1,category1, category1.label,specInfo1=Info1, myylim)
    return(p1)
  })
  names(p1.list) = genes.plot
  if(!is.null(x2)){
    if((!is.data.frame(x2$data))|(!is.data.frame(x2$rhythm))){
      stop("x2: data and rhythmicity parameters must be dataframes")
    }
    if(length(x2$data)!=length(x2$tod)){
      stop("x2: Number of samples in data does not match that in tod. ")
    }
    if(nrow(x2$data)!=length(x2$label)){
      stop("x2: Number of labels does not match number of feature in data. ")
    }
    if(!is.null(category1)){
      if(ncol(x2$data)!=length(category1)){
        stop("length of category 1 does not equal number of samples in x2")
      }
    }
    if(is.null(category2.label)){
      category2.label = "category2"
    }

    p2.list0 =
      lapply(1:length(genes.plot), function(a){
        gene2 = genes.plot[a]
        tod2 = x2$tod
        expr1 = as.numeric(x1$data[match(gene2, x1$label),])
        expr2 = as.numeric(x2$data[match(gene2, x2$label),])
        apar2 = x2$rhythm[match(gene2, x2$label),]
        myylim = c(min(expr1, expr2), max(expr1, expr2))
        p1 = suppressMessages(p1.list[[a]]+ggplot2::ylim(myylim[1], myylim[2]))
        p2 = circadianDrawing_one(tod2, expr2, apar2, gene2,category2, category2.label,specInfo1=Info2, myylim)
        return(list(p1 = p1,
                    p2 = p2))
      })
    p1.list = lapply(p2.list0, `[[`, 1)
    p2.list = lapply(p2.list0, `[[`, 2)
    names(p1.list) = names(p2.list) = genes.plot
  }

  if(!is.null(filename)){
    pdf(filename, width, height)
    if(is.null(x2)){
      lapply(1:length(p1.list), function(a){
        print(p1.list[[a]])
      })
    }else{
      lapply(1:length(p1.list), function(a){
        gridExtra::grid.arrange(p1.list[[a]], p2.list[[a]], ncol = 2)
      })
    }

    dev.off()
  }

  if(is.null(x2)){
    return(p.list = list(x1 = p1.list))
  }else{
    return(p.list = list(x1 = p1.list,
                         x2 = p2.list))
  }

}

# pdf(paste0(plot.dir, "/", target.genes[i], ".pdf"), width = 12)
# tod1 = clinical.l1.g1$Corrected.TOD; tod2 = clinical.l1.g2$Corrected.TOD
# expr1 = as.numeric(data.g1[i, ]); expr2 = as.numeric(data.g2[i, ])
# apar1 = rhythm.g1[i, ]; apar2 = rhythm.g2[i, ]
# changeInfo = paste("p(d_sigma2)=", round(file$p_d_sigma2[i], 2), "\np(d_A.sigma2)=", round(file$p_d_A_sigma2[i], 2))
# circadianDrawing_pair(tod1, expr1, apar1, tod2, expr2, apar2,  changeInfo,
#                       specInfo0=paste0(group1,  "_", comparison.groups$con[a]),  specInfo1=levels.group2[1],  specInfo2=levels.group2[2])
# dev.off()

# tod1 = data1.rhythm$tod; tod2 = data2.rhythm$tod
# expr1 = as.numeric(data1.rhythm$data[3859, ]); expr2 = as.numeric(data2.rhythm$data[3859, ])
# apar1 = data1.rhythm$rhythm[3859, ]; apar2 = data2.rhythm$rhythm[3859, ]
# category1 = clinical$Sex.Label[clinical$Diagnosis.3Grp=="CONTROL"]
# category1.label = "gender"
# myylim = c(5, 10)
# specInfo1="atest"
circadianDrawing_one = function(tod1, expr1, apar1, gene1, category1, category1.label,
                                specInfo1=NULL, myylim){

  geneName <- gene1

  amain1 <- paste0(specInfo1, " " ,geneName,
                   "\n",  " p=", round(apar1$pvalue, 2), ", ", "R2= ", round(apar1$R2, 2), ", A= ", round(apar1$A, 2), ", peak= ", round(apar1$peak, 2),
                   "\n", "sigma2= ", round(apar1$sigma2, 2))


  funx1 = function(x){
    apar1$M+apar1$A*cos(2*pi/apar1$P*x+apar1$phase)
  }
  y1 = funx1(tod1)

  df = data.frame(TOD = tod1, Expression = expr1)
  p1 = ggplot2::ggplot(df, ggplot2::aes(x = TOD, y = Expression))+
    ggplot2::geom_point(ggplot2::aes(color = category1), size = 2)+
    ggplot2::geom_function(fun = funx1, color = "red", size = 2)+
    ggplot2::ylim(min(myylim[1], min(y1)), max(myylim[2], max(y1)))+
    ggplot2::ggtitle(amain1)+
    ggplot2::guides(color=ggplot2::guide_legend(title=category1.label))+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="bottom")

  return(p1)
}

# pdf(NULL)
# dev.control(displaylist="enable")
# plot(tod1,expr1, pch=16,cex=2, col = as.factor(category1),
#      main=amain1,xlim=c(-6,18), #ylim = myylim,
#      xlab='TOD',ylab='Expression')
# #legend("topleft", legend = unique(factor(labels)), col = unique(labelColor), pch = 16)
# smoothingSpline1 = smooth.spline(times, pred1, spar=0.35)
# lines(smoothingSpline1,col='red',lwd=4)
# box(which = "plot", lty = "solid",lwd=3)
# p1.base <- recordPlot()
# invisible(dev.off())
#
# # Display the saved plot
# grid::grid.newpage()
# p1.base

circadianDrawing_pair = function(tod1, expr1, apar1, tod2, expr2, apar2,  changeInfo = NULL,
                                 specInfo0=NULL,  specInfo1=NULL,  specInfo2=NULL){
  getPred <- function(parS, xx) {
    #  parS$A * sin(2*pi/24 * (xx + parS$phase)) + parS$offset
    parS$A * cos(2*pi/24 * xx + parS$phase) + parS$M
  }

  stopifnot("The gene labels from two groups do not match." = apar1$genes==apar2$genes)
  geneName <- apar1$genes

  amain1 <- paste0(specInfo0, " ", geneName,
                   "\n", specInfo1, " p=", round(apar1$pvalue, 2), ", ", "R2= ", round(apar1$R2, 2), ", A= ", round(apar1$A, 2), ", peak= ", round(24-24*apar1$phase/(2*pi)),
                   "\n", "sigma2= ", round(apar1$sigma2, 2), " A.sigma = ", round(apar1$A/sqrt(apar1$sigma2), 2))
  amain2 <- paste0(specInfo0, " ", geneName,
                   "\n", specInfo2, " p=", round(apar2$pvalue, 2), ", ", "R2= ", round(apar2$R2, 2), ", A= ", round(apar2$A, 2), ", peak= ", round(24-24*apar2$phase/(2*pi)),
                   "\n", "sigma2= ", round(apar2$sigma2, 2), " A.sigma = ", round(apar2$A/sqrt(apar2$sigma2), 2))

  times <- seq(-6,18,0.1)
  pred1 <- getPred(apar1,times)
  pred2 <- getPred(apar2,times)

  myylim = c(min(c(expr1, expr2)), max(c(expr1, expr2)))

  par(mfrow = c(1, 2))
  plot(tod1,expr1, pch=16,cex=2, #col = "red",
       main=amain1,xlim=c(-6,18), ylim = myylim,
       xlab='TOD',ylab='Expression')
  #    legend("topleft", legend = unique(factor(labels)), col = unique(labelColor), pch = 16)
  smoothingSpline1 = smooth.spline(times, pred1, spar=0.35)
  lines(smoothingSpline1,col='red',lwd=4)
  text(x = -5, y = myylim[2]- (myylim[2]-myylim[1])*0.05, changeInfo, cex = 0.8, pos = 4)
  box(which = "plot", lty = "solid",lwd=3)

  plot(tod2,expr2, pch=16,cex=2, #col = "red",
       main=amain2,xlim=c(-6,18), ylim = myylim,
       xlab='TOD',ylab='Expression')
  #    legend("topleft", legend = unique(factor(labels)), col = unique(labelColor), pch = 16)
  smoothingSpline2 = smooth.spline(times, pred2, spar=0.35)
  lines(smoothingSpline2,col='red',lwd=4)
  box(which = "plot", lty = "solid",lwd=3)
}
