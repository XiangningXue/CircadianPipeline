#' Title
#'
#' @param x1 From gorup 1: a list containing the following component: data: dataframe with rows of genes and columns of samples; tod: time of death or time of expression corresponding to the columns in data component; label: the gene names or other labels of the genes
#' @param x2 From group 2.
#' @param x.joint A list with two components: Rhythmic.Both contains the labels of genes that are rhythmic in both regions; Rhythmic.GT.One contains the labels of genes that are rhythmic in at least one region
#' @param diffPar.method A string. One of  "circacompare", "permutation".
#' @param diffPar.adjust FDR controlling function, one of adj_minP, adj_BH_separate, adj_BH_pool, adj_AWFisher
#' @param diffR2.method A string. One of "LR", "permutation", "bootstrap"
#' @param nSampling A numeric value. Number of permutation performer, only required when method = "permutation". The smallest possible p-value in the result will be 1/(n_gene*nSampling), but larger nSampling takes longer computing time.
#' @param Sampling.save A local directory where you want to save the permutation result for future use. If NULL then no result will be saved.
#' @param Sampling.file.label pecial label for file name of permutation if you want to save the permutation file
#' @param p.adjust.method choose valid input for p.adjust.method by checking p.adjust.methods for p.adjust() in R package `stat`
#' @param alpha  number between 0 to 1. The critical level used in differential parameter test.
#' @param parallel TRUE/FALSE. If TRUE, parallel computing using mclapply will be used, which does not work on windows system.
#' @param cores number of cores used if using parallel computing
#'
#' @return A list with two components: diffPar: result of differential parameter test for genes rhythmic in both groups; diffR2: result of differential R2 test for genes rhythmic in at least one group
#' @export
#'
#' @examples
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
#' x1 = list(data = as.data.frame(noise.mat1 + signal.mat1),
#' tod = x1.time, label = paste("gene", seq_len(100)))
#' x1.Rhythm = CP_Rhythmicity(x1, parallel = FALSE)
#' x2 = list(data = as.data.frame(noise.mat2 + signal.mat2),
#' tod = x2.time, label = paste("gene", seq_len(100)))
#' x2.Rhythm = CP_Rhythmicity(x2, parallel = FALSE)
#'
#' x.joint = CP_JointRhythmicity(x1.Rhythm, x2.Rhythm, method = "AW-Fisher",
#' para.list = list(param = "pvalue", range = "less", value = 0.05))
#'
#' diffRhythm.res = CP_DiffRhythmicity(x1.Rhythm, x2.Rhythm, x.joint)
#'

CP_DiffRhythmicity = function(x1 = data1.rhythm, x2 = data2.rhythm, x.joint = joint.rhythm,
                              diffPar.method = "circacompare", diffPar.adjust = adj_minP,
                              diffR2.method = "LR",
                              nSampling=1000, #permutation.load = TRUE,
                              Sampling.save = getwd(), Sampling.file.label = "Group1",
                              alpha = 0.05, p.adjust.method = "BH", parallel = FALSE, cores = 5){

  #checking
  if((!is.data.frame(x1$data))|(!is.data.frame(x2$data))){
    stop("x1$data and x2$data must be a dataframe")
  }
  if((ncol(x1$data)!=length(x1$tod))|(ncol(x2$data)!=length(x2$tod))){
    stop("Number of samples in data does not match that in number of tod. ")
  }
  if((nrow(x1$data)!=length(x1$label))|(nrow(x2$data)!=length(x2$label))){
    stop("Number of labels does not match number of genes in data. ")
  }


  all.genes = unique(unlist(x.joint))
  if((!all(all.genes%in%x1$label))|(!all(all.genes%in%x2$label))){
    stop("Genes in x.joint are not contained in both x1 and x2.")
  }

  p1 = x1$rhythm$P[match(all.genes, x1$label)]
  p2 = x2$rhythm$P[match(all.genes, x2$label)]
  period = p1[1]
  if(!(all(c(p1, p2)==period))){
    stop("The period is not the same for all the genes in x1 and x2. ")
  }

  adjust.phase = function(phase = phase, pp = 2*pi){
    phase = phase%%(pp)
    return(phase)
  }


  if(diffPar.method == "permutation"&diffR2.method== "permutation"){
    if((!is.null(Sampling.save))&!dir.exists(Sampling.save)){
      dir.create(file.path(Sampling.save), recursive = TRUE)
      message(paste0("Directory has been created. Permutation results will be saved in ", Sampling.save))
    }
    diff.tab = diff_rhythmicity_permutation(x1, x2, x.joint$Rhythmic.GT.One, period, nSampling, Sampling.save, Sampling.file.label, parallel, cores)
    diffPar.tab = diff.tab[match(x.joint$Rhythmic.Both, diff.tab$label), ]
    diffR2.tab = diff.tab[, c("delta.R2", "p.R2", "label")]
  }else{
    #differential parameter test
    #"two_cosinor1", "two_cosinor2". The method "two_cosinor1" will test the overlap of two CI, and method "two_cosinor2" test if the estimate of group 2 is in CI of group 1.
    overlap.g = x.joint$Rhythmic.Both
    if(diffPar.method == "two_cosinor1"){
      x.list = lapply(1:length(overlap.g), function(a){
        one.gene.data = data.frame(time = c(x1$tod, x2$tod),
                                   measure = c(as.numeric(x1$data[match(overlap.g[a], x1$label), ]),
                                               as.numeric(x2$data[match(overlap.g[a], x2$label), ])),
                                   group = c(rep(0, length(x1$tod)), rep(1, length(x2$tod))))
        return(one.gene.data)
      })

      two_cosinor.res = option_parallel(x.list, function(one.data){
        one.res = two_cosinor_OLS(tod = one.data$time, y = one.data$measure, group = one.data$group)
        return(one.res)
      }, parallel, cores)

      diffPar.tab = do.call(rbind.data.frame, lapply(two_cosinor.res, function(one.res){
        as.data.frame(list(delta.M = one.res$g2$M$est-one.res$g1$M$est,
                           delta.A = one.res$g2$A$est-one.res$g1$A$est,
                           delta.phase = one.res$g2$phase$est-one.res$g1$phase$est,
                           delta.peak = one.res$g2$peak-one.res$g1$peak,
                           p.global = one.res$test$global.pval,
                           M.ind = one.res$test$M.ind[1],
                           A.ind = one.res$test$phase.ind[1],
                           phase.ind = one.res$test$phase.ind[1]))
      }))

      diffPar.tab$q.global = stats::p.adjust(diffPar.tab$p.global, p.adjust.method)

    }else if(diffPar.method == "two_cosinor2"){
      x.list = lapply(1:length(overlap.g), function(a){
        one.gene.data = data.frame(time = c(x1$tod, x2$tod),
                                   measure = c(as.numeric(x1$data[match(overlap.g[a], x1$label), ]),
                                               as.numeric(x2$data[match(overlap.g[a], x2$label), ])),
                                   group = c(rep(0, length(x1$tod)), rep(1, length(x2$tod))))
        return(one.gene.data)
      })

      two_cosinor.res = option_parallel(x.list, function(one.data){
        one.res = two_cosinor_OLS(tod = one.data$time, y = one.data$measure, group = one.data$group)
        return(one.res)
      }, parallel, cores)

      diffPar.tab = do.call(rbind.data.frame, lapply(two_cosinor.res, function(one.res){
        as.data.frame(list(delta.M = one.res$g2$M$est-one.res$g1$M$est,
                           delta.A = one.res$g2$A$est-one.res$g1$A$est,
                           delta.phase = one.res$g2$phase$est-one.res$g1$phase$est,
                           delta.peak = one.res$g2$peak-one.res$g1$peak,
                           p.global = one.res$test$global.pval,
                           M.ind = one.res$test$M.ind[2],
                           A.ind = one.res$test$phase.ind[2],
                           phase.ind = one.res$test$phase.ind[2]))
      }))

      diffPar.tab$q.global = stats::p.adjust(diffPar.tab$p.global, p.adjust.method)

    }else if(diffPar.method == "circacompare"){
      x.list = lapply(1:length(overlap.g), function(a){
        one.gene.data = data.frame(time = c(x1$tod, x2$tod),
                                   measure = c(as.numeric(x1$data[match(overlap.g[a], x1$label), ]),
                                               as.numeric(x2$data[match(overlap.g[a], x2$label), ])),
                                   group = factor(c(rep(1, length(x1$tod)), rep(2, length(x2$tod)))))
        return(one.gene.data)
      })

      circa.res = option_parallel(x.list, function(one.data){
        one.res = circacompare::circacompare(one.data, col_time = "time", col_group = "group", col_outcome = "measure", period,
                                             alpha_threshold = 1)
        return(one.res)
      }, parallel, cores)

      diffPar.tab = do.call(rbind.data.frame, lapply(circa.res, function(one.res){
        as.data.frame(list(delta.M = one.res[[2]][6, 2],
                           p.M = one.res[[2]][7, 2],
                           delta.A = one.res[[2]][10, 2],
                           p.A = one.res[[2]][11, 2],
                           delta.phase = adjust.phase(2*pi-one.res[[2]][14, 2]/period),
                           delta.peak = one.res[[2]][14, 2],
                           p.phase = one.res[[2]][15, 2]))
      }))

      diffPar.tab = cbind.data.frame(label = overlap.g, diffPar.tab)
      diffPar.tab = diffPar.adjust(diffPar.tab,  p.adjust.method, alpha)
      # circa.res.plot = lapply(1:length(circa.res), function(a){
      #   p = circa.res[[a]][[1]]+ggtitle(paste0(overlap.g[a]))
      #   return(p)
      # })
    }else if(diffPar.method == "permutation"){
      if((!is.null(Sampling.save))&!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        message(paste0("Directory has been created. Permutation results will be saved in ", Sampling.save))
      }
      diffPar.tab = diff_rhythmicity_permutation(x1, x2, overlap.g, period, nSampling, Sampling.save, Sampling.file.label, parallel, cores)
      diffPar.tab = diffPar.adjust(diffPar.tab,  p.adjust.method, alpha)
    }



    #differential R2 test
    overlap.g = x.joint$Rhythmic.GT.One
    if(diffR2.method == "LR"){
      res.list = lapply(1:length(overlap.g), function(a){
        one.res = differentialR2::LR_deltaR2(tt1 = x1$tod, yy1 = as.numeric(x1$data[a, ]),
                             tt2 = x2$tod, yy2 = as.numeric(x2$data[a, ]),
                             FN = TRUE)
        one.res.tab = data.frame(delta.R2 = x2$rhythm$R2[a]-x1$rhythm$R2[a], pvalue = one.res)
        return(one.res.tab)
      })
      diffR2.tab = do.call(rbind.data.frame, res.list)
      diffR2.tab = cbind.data.frame(label = overlap.g, diffR2.tab)
    }else if(diffR2.method == "permutation"){
      if((!is.null(Sampling.save))&!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        message(paste0("Directory has been created. Permutation results will be saved in ", Sampling.save))
      }
      diff.tab = diff_rhythmicity_permutation(x1, x2, overlap.g, period, nSampling, Sampling.save, Sampling.file.label, parallel, cores)
      diffR2.tab = diff.tab[, c("delta.R2", "p.R2", "label")]
    }else if(method == "bootstrap"){
      if((!is.null(Sampling.save))&!dir.exists(Sampling.save)){
        dir.create(file.path(Sampling.save), recursive = TRUE)
        message(paste0("Directory has been created. Boostrap results will be saved in ", Sampling.save))
      }
      diff.tab = diff_rhythmicity_bootstrap(x1, x2, overlap.g, period, nSampling, Sampling.save, Sampling.file.label, parallel, cores)
      diffR2.tab = diff.tab[, c("delta.R2", "p.R2", "label")]
    }

    diffR2.tab$qvalue = stats::p.adjust(diffR2.tab$pvalue, p.adjust.method)

  }

  return(list(diffPar = diffPar.tab,
              diffR2 = diffR2.tab))
}
