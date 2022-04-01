#' Title Scatter Plots of Expression vs. Time.
#'
#' @param x CP_Rhythmicity() output
#' @param genes.plot a vector of character strings. Names of genes to be plotted. Names should be same in nomenclature as gname in the data list. If NULL, the top 10 most rhythmic genes from group I will be used
#' @param Info1 character string. Used in the plot title for group I.
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param filename character string of the filename for exporting. If NULL, plots are not saved.
#' @param file.width height of the export plot
#' @param file.height width of the export plot
#'
#' @return A list of ggplot2 plots.
#' @export
#'
#' @examples
#' x = CP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' rhythm.res = CP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' rhythm.plots = CP_ScatterPlot(rhythm.res)
#' #to display plot in Rstudio use CP_PlotDisplay()
#' #display the first five plots
#' CP_PlotDisplay(rhythm.plots, id = 1:5)
CP_ScatterPlot = function(x, genes.plot = NULL,
                   Info1 = "gI", Info2 = "gII",
                   filename = NULL, file.width = 8, file.height = 8){
  if(all(c("x1", "x2", "rhythm.joint")%in%names(x))){
    x1 = x$x1
    x2 = x$x2
  }else{
    x1 = x
    x2 = NULL
  }
  if(is.null(x2)){
    stopifnot("x1$data must be dataframe" = is.data.frame(x1$data))
    stopifnot("x1 does not contain result from CP_Rhythmicity()" = !is.null(x1$rhythm))
    stopifnot("Number of samples in data does not match that in time. " = ncol(x1$data)==length(x1$time))
    stopifnot("Please input the gene labels x1$gname. " = !is.null(x1$gname))
    stopifnot("Number of gnames does not match number of genes in data. " = nrow(x1$data)==length(x1$gname))
  }else{
    stopifnot("x1$data must be dataframe and x2$data must be dataframe. " = (is.data.frame(x1$data)&is.data.frame(x2$data)))
    stopifnot("x1 or x2 do not contain result from CP_Rhythmicity()" = !(is.null(x1$rhythm)|is.null(x2$rhythm)))
    stopifnot("Number of samples in data does not match that in time. " = (ncol(x1$data)==length(x1$time))&(ncol(x2$data)==length(x2$time)))
    stopifnot("Please input the gene labels x1$gname and x2$gname. " = !(is.null(x1$gname)|is.null(x2$gname)))
    stopifnot("Number of gnames does not match number of genes in data. " = (nrow(x1$data)==length(x1$gname))&(nrow(x2$data)==length(x2$gname)))
    gname.overlap = intersect(x1$gname, x2$gname)
    stopifnot("There are no overlapping genes between x1$gname and x2$gname. " = length(gname.overlap)>0)
  }

  if(is.null(genes.plot)){
    genes.plot = x1$gname[order(x1$rhythm$pvalue)][1:ifelse(length(x1$gname)>10, 10, length(x1$gname))]
  }
  P = x1$P
  p1.list = lapply(1:length(genes.plot), function(a){
    gene1 = genes.plot[a]
    t1 = x1$time
    expr1 = as.numeric(x1$data[match(gene1, x1$gname),])
    apar1 = x1$rhythm[match(gene1, x1$rhythm$gname),]

    myylim = c(min(expr1, apar1$M-apar1$A), max(expr1, apar1$M+apar1$A))
    # tod1=t1; period=P
    p1 = circadianDrawing_one(t1, expr1, apar1, gene1, P, specInfo1=Info1, myylim)
    return(p1)
  })

  names(p1.list) = genes.plot
  if(!is.null(x2)){

    p2.list0 =
      lapply(1:length(genes.plot), function(a){
        gene2 = genes.plot[a]
        t2 = x2$time
        expr1 = as.numeric(x1$data[match(gene2, x1$gname),])
        expr2 = as.numeric(x2$data[match(gene2, x2$gname),])
        apar1 = x1$rhythm[match(gene2, x1$rhythm$gname),]
        apar2 = x2$rhythm[match(gene2, x2$rhythm$gname),]
        myylim = c(min(expr1, expr2, apar1$M-apar1$A, apar2$M-apar2$A), max(expr1, expr2, apar1$M+apar1$A, apar2$M+apar2$A))
        p1 = suppressMessages(p1.list[[a]]+ggplot2::ylim(myylim[1], myylim[2]))
        p2 = circadianDrawing_one(t2, expr2, apar2, gene2, P, specInfo1=Info2, myylim)
        return(list(p1 = p1,
                    p2 = p2))
      })
    p1.list = lapply(p2.list0, `[[`, 1)
    p2.list = lapply(p2.list0, `[[`, 2)
    names(p1.list) = names(p2.list) = genes.plot
  }

  if(!is.null(filename)){
    pdf(filename, file.width, file.height)
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
    return(p1.list)
  }else{
    return(list(x1 = p1.list,
                x2 = p2.list))
  }

}

#' Title Displaying CP_ScatterPlot outputs
#'
#' @param x CP_ScatterPlot output
#' @param id integer. Indexes for plots to display.
#'
#' @export
CP_PlotDisplay = function(x = CP_ScatterPlot(x, genes.plot = NULL,
                                      Info1 = "gI", Info2 = "gII",
                                      filename = NULL, height = 8, width = 8), id = 1:2){
  if(length(x)==1){
    return(x[id])
  }else{
    return(lapply(id, function(a.g){
      gridExtra::grid.arrange(x[[1]][[a.g]], x[[2]][[a.g]], ncol = 2)
    }))
  }
}

circadianDrawing_one = function(tod1, expr1, apar1, gene1, period,
                                specInfo1=NULL, myylim){

  geneName <- gene1
  fun.cosinor = function(x){
    apar1$M+apar1$A*cos(2*pi/period*x+apar1$phase)
  }
  amain1 <- paste0(specInfo1, " " ,geneName,
                   "\n",  " p=", round(apar1$pvalue, 4), ", ", "R2= ", round(apar1$R2, 2), ", A= ", round(apar1$A, 2), ", peak= ", round(apar1$peak, 1),
                   "\n", "sigma= ", round(apar1$sigma, 2))

  df = data.frame(Time = tod1, Expression = expr1)
  p1 = ggplot2::ggplot(df, ggplot2::aes(x = Time, y = Expression))+
    ggplot2::geom_point(size = 2)+
    ggplot2::geom_function(fun = fun.cosinor, color = "red", size = 2)+
    ggplot2::ylim(myylim[1], myylim[2])+
    ggplot2::ggtitle(amain1)+
    # ggplot2::guides(color=ggplot2::guide_legend(title=category1.label))+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position="bottom")

  return(p1)
}

#' Title Heatmaps for Rhythmicity Signals
#'
#' @param x CP_Rhythmicity() output
#' @param genes.plot a vector of character strings. Names of genes to be plotted. Names should be same in nomenclature as gname in the data list. If NULL, the top 100 most rhythmic genes from group I will be used
#' @param col_breaks color breaks
#' @param col_low_mid_high color choice for the heatmaps
#' @param Info1 character string. Used in the plot title for group I.
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param ... other argument to pass to \code{ComplexHeatmap::Heatmap()}
#' @param filename character string of the filename for exporting. If NULL, plots are not saved.
#' @param file.width height of the export plot
#' @param file.height width of the export plot
#'
#' @export
#'
#' @examples
#' x = CP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/4, pi/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' rhythm.res = CP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' CP_PlotHeatmap(rhythm.res)
CP_PlotHeatmap = function(x, genes.plot = NULL,
                          col_breaks = c(seq(-2,0,length=100),
                                         seq(.1,1,length=200),
                                         seq(1.1,2, length = 200)),
                          col_low_mid_high = c("blue","yellow","gold"),
                          Info1 = "group I", Info2 = "group II", ...,
                          filename = NULL, file.width = 8, file.height = 4){
  standardize = function(gene.i){
    x = as.numeric(gene.i)
    names(x) = names(gene.i)
    s.x = (x-mean(x))/sd(x)
    s.x[which(s.x>2)] = 2
    s.x[which(s.x<(-2))] = -2
    return(s.x)
  }

  if(all(c("x1", "x2", "rhythm.joint")%in%names(x))){
    x1 = x$x1
    x2 = x$x2
  }else{
    x1 = x
    x2 = NULL
  }

  if(is.null(genes.plot)){
    genes.plot = x1$rhythm$gname[order(x1$rhythm$pvalue)][1:ifelse(nrow(x1$rhythm)>100, 100, nrow(x1$rhythm))]
    warning()
  }
  stopifnot("x1 is not given" = !is.null(x1))

  if(is.null(x2)){
    if(!all(genes.plot%in%x1$gname)){
      warning("Not all input genes are present in x1. The heatmap will not show the missing genes. ")
    }

    genes.plot = genes.plot[genes.plot%in%x1$gname]
    rhythm.x1 = x1$rhythm[match(genes.plot, x1$rhythm$gname), ]
    mat.x1 = x1$data[match(genes.plot, x1$gname), ]
    rownames(mat.x1) = genes.plot
    mat.x1 = mat.x1[order(rhythm.x1$peak), order(x1$time)]
    mat.x1 = t(apply(mat.x1, 1, standardize))

    p1 = ComplexHeatmap::Heatmap(mat.x1, name = Info1, col = circlize::colorRamp2(col_breaks, colorRampPalette(col_low_mid_high)(n = 500)),
                                 cluster_rows = FALSE, cluster_columns = FALSE, ...)
    if(!is.null(filename)){
      pdf(filename, width = file.width, height = file.height)
      print(p1)
      dev.off()
    }

    return(p1)
  }else{
    if(!(all(genes.plot%in%x1$gname)&all(genes.plot%in%x2$gname))){
      warning("Not all input genes are present in x1 or x2. The heatmap will not show genes in neither x1 nor x2. The side by side heatmaps might contain missingness.")
    }
    genes.plot = genes.plot[genes.plot%in%x1$gname|genes.plot%in%x2$gname]
    rhythm.x1 = x1$rhythm[match(genes.plot, x1$rhythm$gname), ]
    mat.x1 = x1$data[match(genes.plot, x1$gname), ]
    rownames(mat.x1) = genes.plot
    mat.x1 = mat.x1[order(rhythm.x1$peak), order(x1$time)]
    mat.x1 = t(apply(mat.x1, 1, standardize))

    rhythm.x2 = x2$rhythm[match(genes.plot, x2$rhythm$gname), ]
    mat.x2 = x2$data[match(genes.plot, x2$gname), ]
    rownames(mat.x2) = genes.plot
    mat.x2 = mat.x2[order(rhythm.x1$peak), order(x2$time)]
    mat.x2 = t(apply(mat.x2, 1, standardize))

    p1 = ComplexHeatmap::Heatmap(mat.x1, name = Info1, col = circlize::colorRamp2(col_breaks, colorRampPalette(col_low_mid_high)(n = 500)),
                                 cluster_rows = FALSE, cluster_columns = FALSE, ...)
    p2 = ComplexHeatmap::Heatmap(mat.x2, name = Info2, col = circlize::colorRamp2(col_breaks, colorRampPalette(col_low_mid_high)(n = 500)),
                                 cluster_rows = FALSE, cluster_columns = FALSE, ...)
    if(!is.null(filename)){
      pdf(filename, width = file.width, height = file.height)
      print(p1+p2)
      dev.off()
    }
    return(p1+p2)
  }


  # grab_grob <- function(){
  #   gridGraphics::grid.echo()
  #   grid::grid.grab()
  # }
  # heatmap.plots = lapply(list(mat.x1, mat.x2), function(a){
  #   gplots::heatmap.2(mat.x1, Rowv = NA, Colv = NA, main = paste0("\n\n", Info1),
  #                     ylab = info.gene,  col = colorRampPalette(c("blue","yellow","gold"))(n = 499), symkey=FALSE, trace = "none",
  #                     cexRow = 0.5, key = TRUE,density.info = "none", breaks = col_breaks, dendrogram = "none", margins = c(6, 7),
  #                     lmat = rbind(c(4, 4, 3, 3),c(2, 1, 1, 1)),lhei = c(1, 4), lwid = c(0.5,0.5, 1,2))
  #   grab_grob()
  # })

}

#' Title Circos Plots for Phase Summary
#'
#' @param x output from CP_DiffPar()
#' @param type character string. Choices are "square", "circos_connected", "circos_diff", "circos_hist", "circos_hist2".
#' @param sig.cut A list. For single group plot, only genes satisfying sig.cut will be plotted; for two-group plots, all RhyBoth genes will be plotted and genes satisfying sig.cut will be plotted in a different color. The list has three components: \itemize{
#' \item param parameter used for the cutoff. Should be a column in x.
#' \item fun character string. Either "<", or ">"
#' \item val numeric. Number used for the cutoff}
#' @param color.df data frame for color coding. Only used for "circos_diff'. Only one of color.df and sig.cut should be specified. The data frame should have two columns named color and label. label will be displayed as legend. The color vector should have the same order as genes in x$diffPar.tab
#' @param Info1 character string. Used in the plot title for group I.
#' @param Info2 character string. Used in the plot title for group II (if exist).
#' @param filename character string of the filename for exporting. If NULL, plots are not saved.
#' @param file.width height of the export plot
#' @param file.height width of the export plot
#' @param phase.start numeric. What time among the periodic clock do you want the phase start? Default is -6, which is midnight if time is in ZT scale.
#' @param concordance.ref numeric. The radius where the concordance reference be plotted away from \eqn{\Delta}peak = 0.
#' @param single.binwidth numeric. The binwidth for plotting peak histogram. Used in "circos_hist" and "square".
#' @param cir.x.breaks numeric. A vector for breaks for the radius (peak difference) used in "circos_diff". Should only contains values from -period/2 to period/2 and it is recommended that the break is equal spaced.
#' @param cir.y.breaks numeric. A vector for breaks for the angles. Should start with phase.start and end with phase.start+period
#' @param axis.text.size numeric. Size for the axis text.
#' @param legend.position "left”,“top”, “right”, “bottom”, or "none"
#'
#' @export
#'
#' @examples
#' x = CP_sim_data(ngene=1000, nsample=30, A1=c(1, 3), A2=c(1, 3),
#' phase1=c(0, pi/4), phase2=c(pi/2, pi*3/2),
#' M1=c(4, 6), M2=c(4, 6), sigma1=1, sigma2=1)
#'
#' rhythm.res = CP_Rhythmicity(x1 = x[[1]], x2 = x[[2]])
#' #notice that for CP_PlotPhase() the parameters tested by CP_DiffPar must contain "phase"
#' rhythm.diffPar = CP_DiffPar(rhythm.res, "phase")
#' #make scatter plots of \eqn{\delta}peak vs. peak in group I
#' CP_PlotPhase(rhythm.diffPar, "circos_diff")
#' #make a connected phase plot
#' CP_PlotPhase(rhythm.diffPar, "circos_connected")
CP_PlotPhase = function(x,  type,
                        sig.cut = list(param = "pvalue", fun = "<", val = "0.05"),
                        color.df = NULL,
                        Info1 = "groupI", Info2 = "groupII",
                        filename = "NULL", file.width = 8, file.height = 8,
                        phase.start = -6, concordance.ref = 4,
                        single.binwidth = 1,
                        # cir.y.grid = seq(-6,18, 4),
                        cir.x.breaks = seq(-12,12, 4), cir.y.breaks = seq(-6,18, 4),
                        axis.text.size = 12, legend.position="right"){

  adjust.circle = function(x, a.min = -6,  period = 24){
    a.max = a.min + period
    x[x<a.min] = x[x<a.min]+period
    x[x>a.max] = x[x>a.max]-period
    return(x)
  }

  # to.df.radar = function(peak.vec = peak.df$peak, bin.width = 1, phase.start){
  #   a.break = seq(phase.start, phase.start+period, by = single.binwidth)
  #   ranges = cut(peak.vec, breaks = a.break)
  #   a.df = as.data.frame(table(ranges))
  #   return(a.df)
  # }

  stopifnot("Input should be output of CP_Rhythmicity (for single group plot) or CP_DiffPar (for phase difference plot)" = !is.null(x$P))
  stopifnot("type should be one of 'square', 'circos_hist', 'circos_hist2', 'circos_connected', 'circos_diff'. " = type %in% c("square", "circos_connected", "circos_diff", "circos_diff2", "circos_hist", "circos_hist2"))
  stopifnot("At least one of sig.cut and color.df should be NULL" = sum(is.null(sig.cut), is.null(color.df))>=1)

  period = x$P

  if(is.null(x$diffPar.tab)){
    if(!is.null(sig.cut)){
      xx = ifelse(sig.cut$param == "delta.peak", abs(x$rhythm[, sig.cut$param]), x$rhythm[, sig.cut$param])
      peak.df = x$rhythm[.Primitive(sig.cut$fun)(xx, sig.cut$val), ]
    }else{
      peak.df = x$rhythm
    }
    peak.df$peak = adjust.circle(peak.df$peak, phase.start, period)
    pp.file = paste0(filename, "_", Info1, "_", type)
    if(type == "circos_hist"|type == "circos_hist2"|type == "circos_diff"|type == "circos_connected"){
      pp = ggplot2::ggplot(data = peak.df,ggplot2::aes(y = peak))+
        ggplot2::geom_histogram(binwidth = single.binwidth, fill = "#3374b0")  +
        # ggplot2::geom_text(data=data.frame(x=cir.x.breaks, y=sum(cir.y.breaks[1:2])/2, label=cir.x.breaks), ggplot2::aes(x=x, y=y, label = label), nudge_x = -0.2, size=axis.text.size*1/3) +
        ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max),
                                    labels = cir.y.breaks) +
        ggplot2::xlab("") + ggplot2::ylab(paste0("Peak time in ", Info1)) +
        ggplot2::ggtitle(paste0("Peak histogram of ", Info1))+
        ggplot2::theme_bw()+
        ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = axis.text.size),
                       axis.text.y = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                       plot.title = ggplot2::element_text(size=axis.text.size+2)) +
        ggplot2::coord_polar(theta="y", start=a.min)

      # df.radar0 = to.df.radar(peak.df$peak, bin.width = 1, phase.start)
      # df.radar = data.frame(max = max(df.radar0$Freq), min = 0, `Info1` = df.radar0$Freq)
      # rownames(df.radar) = df.radar0$ranges
      # df.radar = as.data.frame(t(df.radar))
      #   fmsb::radarchart(df.radar)
    }else if(type == "square"){
      pp = ggplot2::ggplot(data = peak.df,ggplot2::aes(x = peak))+
        ggplot2::geom_histogram(binwidth = single.binwidth, fill = "#3374b0")  +
        # ggplot2::geom_text(data=data.frame(x=cir.x.breaks, y=sum(cir.y.breaks[1:2])/2, label=cir.x.breaks), ggplot2::aes(x=x, y=y, label = label), nudge_x = -0.2, size=axis.text.size*1/3) +
        ggplot2::scale_x_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max),
                                    labels = cir.y.breaks) +
        ggplot2::ylab("count") + ggplot2::xlab(paste0("Peak time in ", Info1)) +
        ggplot2::ggtitle(paste0("Peak histogram of ", Info1))+
        ggplot2::theme_bw()+
        ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = axis.text.size),
                       axis.text.y = ggplot2::element_text(size = axis.text.size),
                       panel.border = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                       plot.title = ggplot2::element_text(size=axis.text.size+2))
    }
  }else{

    # phase.start = -6; concordance.ref = 4;
    # single.binwidth = 1;
    # cir.y.grid = seq(-6,18, 4); cir.x.breaks = seq(-12,12, 4); cir.y.breaks = seq(-6,18, 4);
    # axis.text.size = 12; legend.position="right"

    x = x$diffPar.tab
    a.min = phase.start
    a.max = phase.start+period
    peak1 = x$peak1
    peak2 = x$peak2


    peak.df = data.frame(peak1 = peak1, peak2 = peak2)
    peak.df$peak1 = adjust.circle(peak.df$peak1, phase.start, period)
    peak.df$peak2 = adjust.circle(peak.df$peak2, phase.start, period)
    peak.df$delta.peak = peak.df$peak2 - peak.df$peak1

    peak.df$delta.peak[peak.df$delta.peak>(period/2)] = peak.df$delta.peak[peak.df$delta.peak>(period/2)] -period
    peak.df$delta.peak[peak.df$delta.peak< -(period/2)] = peak.df$delta.peak[peak.df$delta.peak< -(period/2)] +period
    pp.file = paste0(filename, "_", Info1, "_", Info2, "_", type)

    if(!is.null(sig.cut)){
      xx = ifelse(rep(sig.cut$param == "delta.peak", nrow(x)), abs(x[, sig.cut$param]), x[, sig.cut$param])
      sig.color = .Primitive(sig.cut$fun)(xx, sig.cut$val)
      sig.color = ifelse(sig.color, "sig", "none") #red: #b33515
      legend.label = ifelse(sig.cut$param == "delta.peak", paste0("|peak difference|", sig.cut$fun, sig.cut$val), c(paste(unlist(sig.cut), collapse="")))
      #paste(unlist(sig.cut), collapse="")
    }else if(!is.null(color.df)){
      sig.color = color.df$label
      sig.color.breaks = names(table(color.df$label))
      sig.color.values.ind = apply(table(color.df$color, color.df$label), 2, function(a){which(a!=0)})
      if(class(sig.color.values.ind)=="list"){
        stop("Please make sure that label and color in color.df are one-to-one matched. ")
      }
      sig.color.values = rownames(table(color.df$color, color.df$label))[sig.color.values.ind]
    }else{
      sig.color = NULL
    }

    if(type == "circos_diff"){
      # cir.x.breaks = seq(-12,12, 4)
      # cir.x.breaks = c(4, 8, 12, -8, -4, 0, 4)
      peak.df$delta.peak2 = ifelse(peak.df$delta.peak>=cir.x.breaks[1], peak.df$delta.peak, peak.df$delta.peak+period)
      cir.x.breaks2 = ifelse(cir.x.breaks>=cir.x.breaks[1], cir.x.breaks, cir.x.breaks+period)
      if(cir.x.breaks[1]==tail(cir.x.breaks, 1)){cir.x.breaks2[length(cir.x.breaks2)] = cir.x.breaks2[length(cir.x.breaks2)] + period}
      highlight.radius = period/6
      if(sum(cir.x.breaks2%%period==0)){#if 0 is in the given axis
        highlight.center =  cir.x.breaks2[which(cir.x.breaks2%%period==0)]
      }else{
        cir.x.breaks2.resid = cir.x.breaks2-period
        cir.x.breaks2.resid.min.ind = which.min(abs(cir.x.breaks2.resid))
        highlight.center = cir.x.breaks2[cir.x.breaks2.resid.min.ind]-cir.x.breaks2.resid[cir.x.breaks2.resid.min.ind]
      }

      InRange = function(x, y){
        x.min = min(x); x.max = max(x)
        if(y<=x.max&y>=x.min){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }
      if(InRange(cir.x.breaks2, highlight.center+highlight.radius)&InRange(cir.x.breaks2, highlight.center-highlight.radius)){
        rects = data.frame(start = highlight.center-highlight.radius, end = highlight.center+highlight.radius, group = 1)
      }else if(InRange(cir.x.breaks2, highlight.center+highlight.radius)){
        rects1 = data.frame(start = min(cir.x.breaks), end = highlight.center+highlight.radius, group = 1)
        rects2 = data.frame(start = highlight.center-highlight.radius+period, end = max(cir.x.breaks), group = 2)
        rects = rbind.data.frame(rects1, rects2)
      }else{
        InRange(cir.x.breaks2, highlight.center-highlight.radius)
        rects1 = data.frame(start = highlight.center-highlight.radius, end = max(cir.x.breaks), group = 1)
        rects2 = data.frame(start = min(cir.x.breaks), end = highlight.center+highlight.radius-period, group = 2)
        rects = rbind.data.frame(rects1, rects2)
      }

      pp = ggplot2::ggplot(data = peak.df,ggplot2::aes(x = delta.peak2, y = peak1))+
        ggplot2::geom_rect(data=rects, inherit.aes=FALSE, ggplot2::aes(xmin=start, xmax=end, ymin=min(cir.y.breaks),
                                                                       ymax=max(cir.y.breaks), group=group), color="transparent", fill="darkgreen", alpha=0.1)+
        ggplot2::geom_vline(xintercept = c(cir.x.breaks2[1], tail(cir.x.breaks2, 1)), linetype="dashed", color="grey80") +
        # ggplot2::geom_vline(xintercept = c(-1*concordance.ref,concordance.ref), linetype="dashed", color="darkgreen",size=1, alpha = 0.6) +
        ggplot2::geom_vline(xintercept = highlight.center,  color = "blue", alpha = 0.6) +
        # ggplot2::geom_hline(yintercept = cir.y.grid, color="grey80") +
        ggplot2::geom_point(size=0.4, ggplot2::aes(color = sig.color), alpha = 0.8) +
        {if(!is.null(color.df)) ggplot2::scale_color_manual(name = "color",
                                                            breaks = sig.color.breaks,
                                                            values = sig.color.values,
                                                            labels = sig.color.breaks)}+
        {if(!is.null(sig.cut)) ggplot2::scale_color_manual(name = "color",
                                                           breaks = c("sig"),
                                                           values = c("sig" = "#b33515", "none"= "dark grey"),
                                                           labels=c(paste(unlist(sig.cut), collapse="")))}+
        ggplot2::geom_text(data=data.frame(x=cir.x.breaks2, y=sum(cir.y.breaks[1:2])/2, label=cir.x.breaks), ggplot2::aes(x=x, y=y, label = label), nudge_x = -0.2, size=axis.text.size*1/3) +
        ggplot2::xlab("") + ggplot2::ylab(paste0("Angles: peak time in ", Info1, "\n",
                                                 "Radius: ", "peak difference (", Info2, "-", Info1, ")")) +
        ggplot2::ggtitle(paste0("Peak difference between ", Info1, " and ", Info2))+
        ggplot2::scale_x_continuous(breaks = cir.x.breaks2, limits = c(cir.x.breaks2[1]-1.5, tail(cir.x.breaks2, 1)), expand = c(0,0)) +
        ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max),
                                    labels = cir.y.breaks) +
        ggplot2::theme_bw()+
        ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = axis.text.size),
                       axis.text.y = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                       plot.title = ggplot2::element_text(size=axis.text.size+2)) +
        ggplot2::coord_polar(theta="y", start=a.min)
    }else if(type == "circos_connected"){
      peak.df.long = data.frame(peak = c(peak.df$peak1, peak.df$peak2),
                                delta.peak = rep(peak.df$delta.peak, 2),
                                group = rep(c(1, 2), each = nrow(peak.df)),
                                gene =  rep(seq_along(1:nrow(peak.df)), 2))
      pp = ggplot2::ggplot(data = peak.df.long)+
        ggplot2::geom_line(alpha = 0.3, ggplot2::aes(x = group, y = peak, group = gene, color = abs(delta.peak)))+
        ggplot2::scale_color_gradient(name = "Line color=|Peak difference|", low = "#e8f7ff", high = "#03a1fc")+ #very pale blue to not-so-pale blue
        ggplot2::geom_point(size=2, stroke = 0, alpha = 0.6, shape = 21, ggplot2::aes(x = group, y = peak, color = , fill = rep(sig.color, 2))) + #, ggplot2::aes(color = rep(sig.color, 2))
        {if(!is.null(sig.cut)) ggplot2::scale_fill_manual(name = "Point color",
                                                          breaks = c("sig"),
                                                          values = c("sig" = "#03a1fc", "none"= "dark grey"),
                                                          labels=c(paste(unlist(sig.cut), collapse="")))}+
        ggplot2::geom_text(data=data.frame(x=c(1, 2), y=a.min, label=c(Info1, Info2)), ggplot2::aes(x=x, y=y, label = label), nudge_x = -0.2, size=axis.text.size*1/3) +
        ggplot2::xlab(paste0("")) + ggplot2::ylab("") +
        ggplot2::ggtitle(paste0("Connected peak time between ", Info1, " and ", Info2))+
        ggplot2::scale_x_continuous(breaks=seq(0, 2, 1), limits = c(0, 2), expand = c(0,0)) +
        ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max), labels = cir.y.breaks) +
        ggplot2::theme_bw()+
        ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = axis.text.size),
                       axis.text.y = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.position = legend.position,
                       plot.title = ggplot2::element_text(size=axis.text.size+2)) +
        ggplot2::coord_polar(theta="y", start=a.min)
    }else if(type == "circos_hist"|type == "circos_hist2"){
      peak.df.long = data.frame(peak = c(peak.df$peak1, peak.df$peak2),
                                delta.peak = rep(peak.df$delta.peak, 2),
                                group = factor(rep(c(Info1, Info2), each = nrow(peak.df))))
      pp = ggplot2::ggplot(data = peak.df.long, ggplot2::aes(y = peak, fill = group))+
        #binwidth = single.binwidth, fill = "#3374b0"
        ggplot2::geom_histogram(binwidth = single.binwidth, position = 'identity', alpha = 0.6)+
        ggplot2::xlab(paste0("")) + ggplot2::ylab("") +
        ggplot2::ggtitle(paste0("Histogram of peak time of ", Info1, " and ", Info2))+
        ggplot2::scale_y_continuous(breaks = cir.y.breaks, limits = c(a.min,a.max), labels = cir.y.breaks) +
        ggplot2::theme_bw()+
        ggplot2::theme(aspect.ratio = 1, axis.line = ggplot2::element_blank(),
                       axis.ticks = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = axis.text.size),
                       axis.text.y = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.position = legend.position,
                       plot.title = ggplot2::element_text(size=axis.text.size+2)) +
        ggplot2::coord_polar(theta="y", start=a.min)+
        {if(type == "circos_hist2")ggplot2::facet_wrap(~group)}+
        {if(type == "circos_hist2")ggplot2::theme(legend.position="none")}

    } else if(type == "square"){
      peak.df2 = peak.df
      peak.df2$peak2 =  peak.df2$peak1+peak.df2$delta.peak
      pp = ggplot2::ggplot(data = peak.df2, ggplot2::aes(x = peak1, y = peak2, color = sig.color)) +
        ggplot2::geom_point(alpha = 0.8)+
        {if(!is.null(sig.cut)) ggplot2::scale_color_manual(name = "color",
                                                           breaks = c("sig"),
                                                           values = c("sig" = "#b33515", "none"= "black"),
                                                           labels=c(paste(unlist(sig.cut), collapse="")))}+
        ggplot2::geom_abline(intercept=0,slope=1,colour='blue',size=1)+
        ggplot2::geom_abline(intercept=concordance.ref,slope=1,colour='darkgreen',size=1,linetype=2)+
        ggplot2::geom_abline(intercept=-1*concordance.ref,slope=1,colour='darkgreen',size=1,linetype=2)+
        ggplot2::labs(title = paste0("Peak difference btw. ", Info1, " and ", Info2), x = Info1, y = Info2)+
        ggplot2::theme_bw()+
        ggplot2::theme(axis.ticks.x=ggplot2::element_line(colour="black",size=1,linetype=1),
                       axis.text.x=ggplot2::element_text(colour="black",size=axis.text.size,vjust=0),
                       axis.ticks.y=ggplot2::element_line(colour="black",size=1,linetype=1),
                       axis.text.y=ggplot2::element_text(colour="black",size=axis.text.size,hjust=0),
                       legend.title = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.text = ggplot2::element_text(size=axis.text.size*0.8),
                       legend.position = legend.position,
                       plot.title = ggplot2::element_text(size=axis.text.size+2))
    }
  }
  if(!is.null(filename)){
    pdf(paste0(pp.file, ".pdf"), width = file.width, height = file.height)
    print(pp)
    dev.off()
  }else{
    print(pp)
  }
  return(pp)
}


