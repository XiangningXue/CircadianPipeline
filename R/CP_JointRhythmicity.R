#' Title Joint Rhythmicity of genes from two data sets
#'
#' This function takes the output of CP_Rhythmicity from both data set and test for the joint rhythmicity of overlapping genes of the two data sets.
#'
#' @param x1 rhythm object for data1, should be an output of CP_Rhythmicity
#' @param x2 rhythm object for data2, should be an output of CP_Rhythmicity
#' @param method one of "AW-Fisher" or "cutoff"
#' @param para.list For method = "AW-Fisher", para.list = list(param = "p", range = "less", value = "0.05"). For method = "cutoff", an example is para.list = list(param = c("A", "pvalue"), range = c("greater", "less"), value = c(0.25, 0.05)), which means genes with A>0.25 and pvalue<0.05 are rhythmic genes
#'
#' @return A list with two components: Rhythmic.Both contains the labels of genes that are rhythmic in both regions; Rhythmic.GT.One contains the labels of genes that are rhythmic in at least one region
#' @export
#'
#' @examples
#'
#' #simulate two data sets
#' x1.time = runif(20, min = 0, max = 24)
#' m1 = rnorm(100, 5); A1 = rnorm(100, 3); phase1 = runif(100, min = 0, max = 2*pi); sigma = 1
#' noise.mat1 = matrix(rnorm(100*20, 0, sigma), ncol = 20, nrow = 100)
#' signal.mat1 = t(sapply(1:100, function(a){m1[a]+A1[a]*cos(2*pi/24*x1.time+phase1[a])}))
#'
#' x2.time = runif(20, min = 0, max = 24)
#' m2 = m1; A2 = 1.5*A1; phase2 = runif(100, min = 0, max = 2*pi); sigma = 1
#' noise.mat2 = matrix(rnorm(100*20, 0, sigma), ncol = 20, nrow = 100)
#' signal.mat2 = t(sapply(1:100, function(a){m2[a]+A2[a]*cos(2*pi/24*x2.time+phase2[a])}))
#'
#' x1 = list(data = noise.mat1 + signal.mat1, tod = x1.time, label = paste("gene", seq_len(100)))
#' x1.Rhythm = CP_Rhythmicity(x1, parallel = FALSE)
#' x2 = list(data = noise.mat2 + signal.mat2, tod = x2.time, label = paste("gene", seq_len(100)))
#' x2.Rhythm = CP_Rhythmicity(x2, parallel = FALSE)
#'
#' x.joint = CP_JointRhythmicity(x1.Rhythm, x2.Rhythm, method = "AW-Fisher", para.list = list(param = "pvalue", range = "less", value = 0.05))
#'
CP_JointRhythmicity = function(x1 = data1.rhythm, x2 = data2.rhythm, method = "AW-Fisher", para.list = list(param = "pvalue", range = "less", value = 0.05)){
  #This function takes the rhythmicity results for two data sets and return genes that rhythmic in at least one/two data sets
  #method: one of "AW-Fisher" or "cutoff"
  #para.list: for method = "AW-Fisher", para.list = list(param = "p", range = "less", value = "0.05")
  # for method = "cutoff",   para.list = list(param = c("A", "pvalue"),
  #  range = c("greater", "less"),
  #  value = c(0.25, 0.05))
  # this means genes with A>0.25 and pvalue<0.05 are rhythmic genes
  # para.list = list(param = c("A", "pvalue"),
  #                    range = c("greater", "less"),
  #                  value = c(0.25, 0.05))
  thresholding = function(xx = list(para.tab, range, value)){
    para.tab = xx$para.tab
    range = xx$range
    value = xx$value
    if(range == "greater"){
      x.select = para.tab>value
    }else if(range == "less"){
      x.select = para.tab<value
    }
    return(x.select)
  }

  out = list(Rhythmic.Both = character(),
             Rhythmic.GT.One = character())
  overlap.g = intersect(x1$label, x2$label)

  if(length(overlap.g)==0){
    stop("There is no overlapping genes between the two data sets. Please check labels. ")
  }
  if(method == "AW-Fisher"){
    data1.pvalue = x1$rhythm$pvalue[match(overlap.g, x1$label)]
    data2.pvalue = x2$rhythm$pvalue[match(overlap.g, x2$label)]

    pmatrix = cbind(data1.pvalue,data2.pvalue)

    res = AWFisher::AWFisher_pvalue(pmatrix)
    qvalue <- stats::p.adjust(res$pvalue, "BH")
    meta_p_q = data.frame(pvalue = res$pvalues, qvalue = qvalue, label = overlap.g)

    #now output gene categories
    if(para.list$range=="less"){
      out$Rhythmic.Both = meta_p_q$label[(meta_p_q[, para.list$param]<para.list$value)&(res$weights[,1]==1)&(res$weights[,2]==1)]
      out$Rhythmic.GT.One = meta_p_q$label[(meta_p_q[, para.list$param]<para.list$value)]
    }else if(para.list$range==""){
      out$Rhythmic.Both = meta_p_q$label[(meta_p_q[, para.list$param]>para.list$value)&(res$weights[,1]==1)&(res$weights[,2]==1)]
      out$Rhythmic.GT.One = meta_p_q$label[(meta_p_q[, para.list$param]>para.list$value)]
    }
  }else if(method == "cutoff"){
    paras.tabs = lapply(1:length(para.list$param), function(a){
      return(list(para.tab = data.frame(data1.param = x1$rhythm[match(overlap.g, x1$label), para.list$param[a]],
                                        data2.param = x2$rhythm[match(overlap.g, x2$label), para.list$param[a]]),
                  range = para.list$range[a],
                  value = para.list$value[a]))
    })
    x.select = lapply(paras.tabs, thresholding)
    x.select.data1 = apply(do.call(cbind, lapply(x.select, function(a){a[, 1]})), 1,
                           function(a){purrr::reduce(a, `&`)})
    x.select.data2 = apply(do.call(cbind, lapply(x.select, function(a){a[, 2]})), 1,
                           function(a){purrr::reduce(a, `&`)})
    out$Rhythmic.Both = overlap.g[x.select.data1&x.select.data2]
    out$Rhythmic.GT.One = overlap.g[x.select.data1|x.select.data2]
  }
  return(out)
}
