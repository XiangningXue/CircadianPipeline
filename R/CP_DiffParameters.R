CP_DiffRhythmicity = function(x1 = data1.rhythm, x2 = data2.rhythm, x.joint = joint.rhythm, period=24, method = "circacompare",
                             nPermutation=1000, permutation.load = TRUE,  permutation.save = getwd(), permutation.file.label = "Group1",
                             alpha = 0.05,  p.adjust.method = "BH", parallel = FALSE, cores = 5){
  overlap.g = x.joint$Rhythmic.Both

  if(method == "circacompare"){
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

    res.tab = do.call(rbind.data.frame, lapply(circa.res, function(one.res){
      as.data.frame(list(delta.mesor = one.res[[2]][6, 2],
                         p.mesor = one.res[[2]][7, 2],
                         delta.A = one.res[[2]][10, 2],
                         p.A = one.res[[2]][11, 2],
                         delta.phase = 2*pi-one.res[[2]][14, 2]/period,
                         p.phase = one.res[[2]][15, 2]))
    }))

    # circa.res.plot = lapply(1:length(circa.res), function(a){
    #   p = circa.res[[a]][[1]]+ggtitle(paste0(overlap.g[a]))
    #   return(p)
    # })
  }else if(method == "permutation"){
    if((permutation.save!="NULL")&!dir.exists(permutation.save)){
      dir.create(file.path(permutation.save), recursive = TRUE)
      message(paste0("Directory has been created. Permutation results will be saved in ", permutation.save))
    }
    diff_rhythmicity_permutation(x1, x2, x.joint, period, nPermutation, permutation.save, permutation.file.label, parallel, cores)
  }

  #do correction
  res.tab$p.combined = pmin(3*pmin(res.tab$p.mesor, res.tab$p.A, res.tab$p.phase), 1)
  res.tab$q.combined = p.adjust(res.tab$p.combined, p.adjust.method)
  the.BH.p.cutoff = max(res.tab$p.combined[res.tab$q.combined<alpha])
  res.tab$ind.global = res.tab$q.combined <alpha
  res.tab$ind.mesor = res.tab$p.mesor<the.BH.p.cutoff/3
  res.tab$ind.A = res.tab$p.A<the.BH.p.cutoff/3
  res.tab$ind.phase = res.tab$p.phase<the.BH.p.cutoff/3

  res.tab$label = overlap.g
  return(res.tab)
}
