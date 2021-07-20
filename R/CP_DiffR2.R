CP_DiffR2 = function(x1 = data1.rhythm, x2 = data2.rhythm, x.joint = joint.rhythm, period=24, method = "LR",
                     nSampling=1000, Sampling.load = FALSE,  Sampling.save = "NULL", Sampling.file.label = "Group1",
                     alpha = 0.05,  p.adjust.method = "BH", parallel = FALSE, cores = 5){


  adjust.phase = function(phase = phase, pp = 2*pi){
    phase = phase%%(pp)
    return(phase)
  }

  #subset the target gene set
  overlap.g = x.joint$Rhythmic.GT.One

  if(method == "permutation"){
    if((Sampling.save!="NULL")&!dir.exists(Sampling.save)){
      dir.create(file.path(Sampling.save), recursive = TRUE)
      message(paste0("Directory has been created. Sampling results will be saved in ", Sampling.save))
    }
    diff.tab = diff_rhythmicity_permutation(x1, x2, overlap.g, period, nSampling, Sampling.save, Sampling.file.label, parallel, cores)
    diffR2.tab = diff.tab[, c("delta.R2", "p.R2", "label")]
  }else if(method == "LR"){
    res.list = lapply(1:length(overlap.g), function(a){
      one.res = differentialR2::LR_deltaR2(tt1 = x1$tod, yy1 = as.numeric(x1$data[a, ]),
                                           tt2 = x2$tod, yy2 = as.numeric(x2$data[a, ]),
                                           FN = TRUE)
      one.res.tab = data.frame(delta.R2 = x2$rhythm$R2[a]-x1$rhythm$R2[a], pvalue = one.res)
      colnames(one.res.tab)[colnames(one.res.tab)=="pvalue"]="p.R2"
      return(one.res.tab)
    })
    diffR2.tab = do.call(rbind.data.frame, res.list)
    diffR2.tab = cbind.data.frame(label = overlap.g, diffR2.tab)
  }else if(method == "bootstrap"){
    if((Sampling.save!="NULL")&!dir.exists(Sampling.save)){
      dir.create(file.path(Sampling.save), recursive = TRUE)
      message(paste0("Directory has been created. Sampling results will be saved in ", Sampling.save))
    }
    diff.tab = diff_rhythmicity_bootstrap(x1, x2, overlap.g, period, nSampling, Sampling.save, Sampling.file.label, parallel, cores)
    diffR2.tab = diff.tab[, c("delta.R2", "p.R2", "label")]

  }
  diffR2.tab$q.R2 = stats::p.adjust(diffR2.tab$p.R2, p.adjust.method)
  return(diffR2 = diffR2.tab)

}
