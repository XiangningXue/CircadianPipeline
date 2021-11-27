#' Title
#'
#' @param x1 From gorup 1: a list containing the following component: data: dataframe with rows of genes and columns of samples; tod: time of death or time of expression corresponding to the columns in data component; label: the gene names or other labels of the genes
#' @param x2 From group 2.
#' @param x.joint A list with two components: Rhythmic.Both contains the labels of genes that are rhythmic in both regions; Rhythmic.GT.One contains the labels of genes that are rhythmic in at least one region
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
                              diffPar = "A&phase&M",
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


# diffPar -----------------------------------------------------------------
  overlap.g = x.joint$Rhythmic.Both
  x.list = lapply(1:length(overlap.g), function(a){
    list(x1.time = x1$tod,
         x2.time = x2$tod,
         y1 = as.numeric(x1$data[match(overlap.g[a], x1$label), ]),
         y2 = as.numeric(x2$data[match(overlap.g[a], x2$label), ]))
  })

  if(diffPar=="A&phase&M"){
    test_diffPar = parallel::mclapply(1:length(x.list), function(a){
      test.overall = two_cosinor_OLS_overall(c(x.list[[a]]$x1.time, x.list[[a]]$x2.time),
                                             c(x.list[[a]]$y1, x.list[[a]]$y2),
                                             c(rep(0, length(x.list[[a]]$x1.time)), rep(1, length(x.list[[a]]$x2.time))),
                                             test = "A&phase&M", CI = FALSE)
      testA = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="amplitude")
      testphase = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="phase")
      testM = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="basal")
      one.row = data.frame(label = overlap.g[a],
                      delta.A = testA$amp_2-testA$amp_1,
                      delta.phase = adjust.phase(testphase$phase_2-testphase$phase_1),
                      delta.M = testM$offset_2-testM$offset_1,
                      p.overall = test.overall,
                      post.hoc.A = test.overall<alpha&testA$pvalue<PostHocP(3, alpha, method = "Sidak"),
                      post.hoc.phase = test.overall<alpha&testphase$pvalue<PostHocP(3, alpha, method = "Sidak"),
                      post.hoc.M = test.overall<alpha&testM$pvalue<PostHocP(3, alpha, method = "Sidak")
                      )
      return(one.row)
    }, mc.cores = cores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$q.overall = stats::p.adjust(diffPar.tab$p.overall, p.adjust.method)
  }else if(diffPar=="A&phase"){
    test_diffPar = parallel::mclapply(1:length(x.list), function(a){
      test.overall = two_cosinor_OLS_overall(c(x.list[[a]]$x1.time, x.list[[a]]$x2.time),
                                             c(x.list[[a]]$y1, x.list[[a]]$y2),
                                             c(rep(0, length(x.list[[a]]$x1.time)), rep(1, length(x.list[[a]]$x2.time))),
                                             test = "A&phase", CI = FALSE)
      testA = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="amplitude")
      testphase = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="phase")
      one.row = data.frame(label = overlap.g[a],
                           delta.A = testA$amp_2-testA$amp_1,
                           delta.phase = adjust.phase(testphase$phase_2-testphase$phase_1),
                           p.overall = test.overall,
                           post.hoc.A = test.overall<alpha&testA$pvalue<PostHocP(2, alpha, method = "Sidak"),
                           post.hoc.phase = test.overall<alpha&testphase$pvalue<PostHocP(2, alpha, method = "Sidak")
      )
      return(one.row)
    }, mc.cores = cores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$q.overall = stats::p.adjust(diffPar.tab$p.overall, p.adjust.method)
  }else if(diffPar=="A"){
    test_diffPar = parallel::mclapply(1:length(x.list), function(a){
      testA = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="amplitude")
      one.row = data.frame(label = overlap.g[a],
                           delta.A = testA$amp_2-testA$amp_1,
                           p.A = testA$pvalue
      )
      return(one.row)
    }, mc.cores = cores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$q.A = stats::p.adjust(diffPar.tab$p.A, p.adjust.method)
  }else if(diffPar=="phase"){
    test_diffPar = parallel::mclapply(1:length(x.list), function(a){
      testphase = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="phase")
      one.row = data.frame(label = overlap.g[a],
                           delta.phase = adjust.phase(testphase$phase_2-testphase$phase_1),
                           p.phase = testphase$pvalue
      )
      return(one.row)
    }, mc.cores = cores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$q.phase = stats::p.adjust(diffPar.tab$p.phase, p.adjust.method)
  }else if(diffPar=="M"){
    test_diffPar = parallel::mclapply(1:length(x.list), function(a){
      testM = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="basal")
      one.row = data.frame(label = overlap.g[a],
                           delta.M = testM$offset_2-testM$offset_1,
                           p.M = testM$pvalue
      )
      return(one.row)
    }, mc.cores = cores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$q.M = stats::p.adjust(diffPar.tab$p.M, p.adjust.method)
  }else if(diffPar=="All"){
    test_diffPar = parallel::mclapply(1:length(x.list), function(a){
      test.overall = two_cosinor_OLS_overall(c(x.list[[a]]$x1.time, x.list[[a]]$x2.time),
                                             c(x.list[[a]]$y1, x.list[[a]]$y2),
                                             c(rep(0, length(x.list[[a]]$x1.time)), rep(1, length(x.list[[a]]$x2.time))),
                                             test = "A&phase&M", CI = FALSE)
      test.overall2 = two_cosinor_OLS_overall(c(x.list[[a]]$x1.time, x.list[[a]]$x2.time),
                                             c(x.list[[a]]$y1, x.list[[a]]$y2),
                                             c(rep(0, length(x.list[[a]]$x1.time)), rep(1, length(x.list[[a]]$x2.time))),
                                             test = "A&phase", CI = FALSE)
      testA = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="amplitude")
      testphase = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="phase")
      testM = diffCircadian::LR_diff(x.list[[a]]$x1.time, x.list[[a]]$y1, x.list[[a]]$x2.time, x.list[[a]]$y2, period = 24, FN = TRUE, type="basal")
      one.row = data.frame(label = overlap.g[a],
                           delta.A = testA$amp_2-testA$amp_1,
                           delta.phase = adjust.phase(testphase$phase_2-testphase$phase_1),
                           delta.M = testM$offset_2-testM$offset_1,
                           p.overall = test.overall,
                           p.overall2 = test.overall2,
                           p.A = testA$pvalue,
                           p.phase = testphase$pvalue,
                           p.M = testM$pvalue,
                           post.hoc.A = test.overall<alpha&testA$pvalue<PostHocP(3, alpha, method = "Sidak"),
                           post.hoc.phase = test.overall<alpha&testphase$pvalue<PostHocP(3, alpha, method = "Sidak"),
                           post.hoc.M = test.overall<alpha&testM$pvalue<PostHocP(3, alpha, method = "Sidak"),
                           post.hoc2.A = test.overall2<alpha&testA$pvalue<PostHocP(2, alpha, method = "Sidak"),
                           post.hoc2.phase = test.overall2<alpha&testphase$pvalue<PostHocP(2, alpha, method = "Sidak")

      )
      return(one.row)
    }, mc.cores = cores)
    diffPar.tab = do.call(rbind.data.frame, test_diffPar)
    diffPar.tab$q.overall = stats::p.adjust(diffPar.tab$p.overall, p.adjust.method)
    diffPar.tab$q.overall2 = stats::p.adjust(diffPar.tab$p.overall2, p.adjust.method)
    diffPar.tab$q.A = stats::p.adjust(diffPar.tab$p.A, p.adjust.method)
    diffPar.tab$q.phase = stats::p.adjust(diffPar.tab$p.phase, p.adjust.method)
    diffPar.tab$q.M = stats::p.adjust(diffPar.tab$p.M, p.adjust.method)
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

  return(list(diffPar = diffPar.tab,
              diffR2 = diffR2.tab))
}

PostHocP = function(p = 4, alpha = 0.05, n, method = "Bonferroni"){
  if(method=="Bonferroni"){
    adjusted.p = alpha/p
  }else if(method == "Sidak"){
    adjusted.p = 1-(1-alpha)^(1/p)
  }else{
    stop("Method must be Bonferroni or Sidak")
  }
  return(adjusted.p)
}

