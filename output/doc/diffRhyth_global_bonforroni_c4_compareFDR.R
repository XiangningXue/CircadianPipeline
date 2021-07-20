##differential rhythmicity: global test
#
#speed3
# Set directory and load libraries-----------------------------------------------------------
rm(list = ls())
dir = "/home/xix66/circadian/ThePipeline/diffRhyth_global_bonferroni_noR2change"
setwd(dir)

# Run with the bonferroni->BH->single test --------------------------------
library(parallel)
library(circacompare)

ind.list = list(
  ind.nochange = 1:650,
  ind.change.offset = 651:700,
  ind.change.A = 701:750,
  ind.change.phase = 751:800,
  ind.change.A.phase = 801:850,
  ind.change.A.offset = 851:900,
  ind.change.phase.offset = 901:950,
  ind.change.A.phase.offset = 951:1000
)
ind.change.global = 651:1000
ind.change.offset = unname(do.call(c, ind.list[grepl("offset", names(ind.list))]))
ind.change.A = unname(do.call(c, ind.list[grepl("A", names(ind.list))]))
ind.change.phase = unname(do.call(c, ind.list[grepl("phase", names(ind.list))]))
#simulate data: each data will have the same set of the genes, just repeat the following for 10 times
method.dir = list(m1 = list(name = "1Bonferroni2BH3Individual", 
                            dir = "res_circacompare_1Bonferroni2BH3Individual"), 
                  m2 = list(name = "1BH2Bonferroni", 
                            dir = "res_circacompare_1BH2Bonferroni0"), 
                  m3 = list(name = "BH", 
                            dir = "res_circacompare_BH0"), 
                  m4 = list(name = "AWFisher",
                            dir = "res_circacompare_AWFisher"), 
                  m2b = list(name = "1BH2Bonferroni_divide3", 
                             dir = "res_circacompare_1BH2Bonferroni"), 
                  m3b = list(name = "BH_divide3", 
                             dir = "res_circacompare_BH"))
all.res = mclapply(1:100, function(r){

  one.res = do.call(rbind.data.frame, lapply(method.dir, function(m){
    one.method = readRDS(paste0(dir, "/", m$dir, "/res", "_n", r, ".rds"))
    one.row = data.frame(r = 1,
               method = m$name, 
               FDR.global = sum(one.method$ind.global[-ind.change.global])/sum(one.method$ind.global), 
               FDR.offset = sum(one.method$ind.Mesor[-ind.change.offset])/sum(one.method$ind.Mesor), 
               FDR.A = sum(one.method$ind.amp[-ind.change.A])/sum(one.method$ind.amp),
               FDR.phase = sum(one.method$ind.phase[-ind.change.phase])/sum(one.method$ind.phase), 
               Recall.global = mean(one.method$ind.global[ind.change.global]), 
               Recall.offset = mean(one.method$ind.Mesor[ind.change.offset]), 
               Recall.A = mean(one.method$ind.amp[ind.change.A]), 
               Recall.phase = mean(one.method$ind.phase[ind.change.phase]))
    one.row$F.global = 2*(1-one.row$FDR.global)*one.row$Recall.global/((1-one.row$FDR.global)+one.row$Recall.global)
    one.row$F.offset = 2*(1-one.row$FDR.offset)*one.row$Recall.offset/((1-one.row$FDR.offset)+one.row$Recall.offset)
    one.row$F.A = 2*(1-one.row$FDR.A)*one.row$Recall.A/((1-one.row$FDR.A)+one.row$Recall.A)
    one.row$F.phase = 2*(1-one.row$FDR.phase)*one.row$Recall.phase/((1-one.row$FDR.phase)+one.row$Recall.phase)
    return(one.row)
  }))
  
}, mc.cores = 20)

library(tidyr)
res.tab = do.call(rbind.data.frame, all.res) 
res.tab.FDR = res.tab[, c("r", "method", colnames(res.tab)[grepl("FDR", colnames(res.tab))])]
res.tab.FDR = gather(res.tab.FDR, "test", "FDR", -r, -method)
res.tab.FDR$test = factor(res.tab.FDR$test, levels = c("FDR.global", "FDR.offset", "FDR.A", "FDR.phase"))
res.tab.Recall = res.tab[, c("r", "method", colnames(res.tab)[grepl("Recall", colnames(res.tab))])]
res.tab.Recall = gather(res.tab.Recall, "test", "Recall", -r, -method)
res.tab.Recall$test = factor(res.tab.Recall$test, levels = c("Recall.global", "Recall.offset", "Recall.A", "Recall.phase"))
res.tab.F = res.tab[, c("r", "method", colnames(res.tab)[grepl("F\\.", colnames(res.tab))])]
res.tab.F = gather(res.tab.F, "test", "F", -r, -method)
res.tab.F$test = factor(res.tab.F$test, levels = c("F.global", "F.offset", "F.A", "F.phase"))

library(ggplot2)
p1 = ggplot(data = res.tab.FDR, aes(x = test, y = FDR, fill = method))+
  geom_boxplot() + scale_y_continuous(breaks=seq(0,1,.05))
p2 = ggplot(data = res.tab.Recall, aes(x = test, y = Recall, fill = method))+
  geom_boxplot()+ scale_y_continuous(breaks=seq(0,1,.05))
p3 = ggplot(data = res.tab.F, aes(x = test, y = F, fill = method))+
  geom_boxplot()+ scale_y_continuous(breaks=seq(0,1,.05))
library(gridExtra)
pdf(paste0(dir, "/FDR control compare.pdf"), width = 20)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()


sum(!(is.na(res.tab$FDR.offset)))#362
sum(is.na(res.tab$FDR.offset))#238
table(is.na(res.tab$FDR.offset), res.tab$method)
posum(!(is.na(res.tab$F.A)))#596
sum(is.na(res.tab$F.A))#4
sum(!(is.na(res.tab$F.phase)))#600
sum(is.na(res.tab$F.phase))#0
ind.list = list(
  ind.nochange = 1:650,
  ind.change.offset = 651:700,
  ind.change.A = 701:750,
  ind.change.phase = 751:800,
  ind.change.A.phase = 801:850,
  ind.change.A.offset = 851:900,
  ind.change.phase.offset = 901:950,
  ind.change.A.phase.offset = 951:1000
)
#test for different categories
# Test for different categories of A -------------------------------------------
subcat.res = lapply(1:100, function(r){
  type.vec = do.call(c, lapply(1:length(ind.list), function(a){
    rep(gsub("ind\\.()", "\\1", names(ind.list)[a]), length(ind.list[[a]]))
  }))
  type.level = gsub("ind\\.()", "\\1", names(ind.list))
  
  one.res = do.call(rbind.data.frame, lapply(method.dir, function(m){
    one.method = readRDS(paste0(dir, "/", m$dir, "/res", "_n", r, ".rds"))
    one.sum = do.call(rbind.data.frame, list(
      data.frame(r = 1, method = m$name, type = "offset", 
                 FDR = sum(((one.method$ind.Mesor)&(!one.method$ind.phase)&(!one.method$ind.amp))[-ind.list$ind.change.offset])/sum((one.method$ind.Mesor)&(!one.method$ind.phase)&(!one.method$ind.amp)), 
                 Recall = mean((one.method$ind.Mesor&(!one.method$ind.phase)&(!one.method$ind.amp))[ind.list$ind.change.offset])), 
      data.frame(r = 1, method = m$name, type = "A", 
                 FDR = sum(((!one.method$ind.Mesor)&(!one.method$ind.phase)&(one.method$ind.amp))[-ind.list$ind.change.A])/sum((!one.method$ind.Mesor)&(!one.method$ind.phase)&(one.method$ind.amp)), 
                 Recall = mean(((!one.method$ind.Mesor)&(!one.method$ind.phase)&(one.method$ind.amp))[ind.list$ind.change.A])), 
      data.frame(r = 1, method = m$name, type = "phase", 
                 FDR = sum(((!one.method$ind.Mesor)&(one.method$ind.phase)&(!one.method$ind.amp))[-ind.list$ind.change.phase])/sum((!one.method$ind.Mesor)&(one.method$ind.phase)&(!one.method$ind.amp)), 
                 Recall = mean(((!one.method$ind.Mesor)&(one.method$ind.phase)&(!one.method$ind.amp))[ind.list$ind.change.phase])), 
      data.frame(r = 1, method = m$name, type = "A.phase", 
                 FDR = sum(((!one.method$ind.Mesor)&(one.method$ind.phase)&(one.method$ind.amp))[-ind.list$ind.change.A.phase])/sum((!one.method$ind.Mesor)&(one.method$ind.phase)&(one.method$ind.amp)), 
                 Recall = mean(((!one.method$ind.Mesor)&(one.method$ind.phase)&(one.method$ind.amp))[ind.list$ind.change.A.phase])), 
      data.frame(r = 1, method = m$name, type = "A.offset", 
                 FDR = sum(((one.method$ind.Mesor)&(!one.method$ind.phase)&(one.method$ind.amp))[-ind.list$ind.change.A.offset])/sum((one.method$ind.Mesor)&(!one.method$ind.phase)&(one.method$ind.amp)), 
                 Recall = mean(((one.method$ind.Mesor)&(!one.method$ind.phase)&(one.method$ind.amp))[ind.list$ind.change.A.offset])), 
      data.frame(r = 1, method = m$name, type = "phase.offset", 
                 FDR = sum(((one.method$ind.Mesor)&(one.method$ind.phase)&(!one.method$ind.amp))[-ind.list$ind.change.phase.offset])/sum((one.method$ind.Mesor)&(one.method$ind.phase)&(!one.method$ind.amp)), 
                 Recall = mean(((one.method$ind.Mesor)&(one.method$ind.phase)&(!one.method$ind.amp))[ind.list$ind.change.phase.offset])), 
      data.frame(r = 1, method = m$name, type = "A.phase.offset", 
                 FDR = sum(((one.method$ind.Mesor)&(one.method$ind.phase)&(one.method$ind.amp))[-ind.list$ind.change.A.phase.offset])/sum((one.method$ind.Mesor)&(one.method$ind.phase)&(one.method$ind.amp)), 
                 Recall = mean(((one.method$ind.Mesor)&(one.method$ind.phase)&(one.method$ind.amp))[ind.list$ind.change.A.phase.offset]))
    ))
    one.sum$F.score = 2*(1-one.sum$FDR)*one.sum$Recall/((1-one.sum$FDR)+one.sum$Recall)
    return(one.sum)
  }))
  return(one.res)
})
sum.subcat.tab = do.call(rbind.data.frame, subcat.res)


library(ggplot2)
p1 = ggplot(data = sum.subcat.tab, aes(x = type, y = FDR, fill = method))+
  geom_boxplot() + scale_y_continuous(breaks=seq(0,1,.05))
p2 = ggplot(data = sum.subcat.tab, aes(x = type, y = Recall, fill = method))+
  geom_boxplot()+ scale_y_continuous(breaks=seq(0,1,.05))
p3 = ggplot(data = sum.subcat.tab, aes(x = type, y = F.score, fill = method))+
  geom_boxplot()+ scale_y_continuous(breaks=seq(0,1,.05))
library(gridExtra)
pdf(paste0(dir, "/FDR control compare_subcategories.pdf"), width = 25)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()


# The components of FDR ---------------------------------------------------
#what are the components of FDR? 
#For none genes, could be either categories of changes
#For offset genes, could be all other categories of changes. 
ind.list = list(
  ind.nochange = 1:650,
  ind.change.offset = 651:700,
  ind.change.A = 701:750,
  ind.change.phase = 751:800,
  ind.change.A.phase = 801:850,
  ind.change.A.offset = 851:900,
  ind.change.phase.offset = 901:950,
  ind.change.A.phase.offset = 951:1000
)
ind.change.global = 651:1000
ind.change.offset = unname(do.call(c, ind.list[grepl("offset", names(ind.list))]))
ind.change.A = unname(do.call(c, ind.list[grepl("A", names(ind.list))]))
ind.change.phase = unname(do.call(c, ind.list[grepl("phase", names(ind.list))]))
#simulate data: each data will have the same set of the genes, just repeat the following for 10 times
method.dir = list(m1 = list(name = "1Bonferroni2BH3Individual", 
                            dir = "res_circacompare_1Bonferroni2BH3Individual"), 
                  m2 = list(name = "1BH2Bonferroni", 
                            dir = "res_circacompare_1BH2Bonferroni0"), 
                  m3 = list(name = "BH", 
                            dir = "res_circacompare_BH0"), 
                  m4 = list(name = "AWFisher",
                            dir = "res_circacompare_AWFisher"), 
                  m2b = list(name = "1BH2Bonferroni_divide3", 
                             dir = "res_circacompare_1BH2Bonferroni"), 
                  m3b = list(name = "BH_divide3", 
                             dir = "res_circacompare_BH"))
type.vec = do.call(c, lapply(1:length(ind.list), function(a){
  rep(gsub("ind\\.()", "\\1", names(ind.list)[a]), length(ind.list[[a]]))
}))
type.level = gsub("ind\\.()", "\\1", names(ind.list))

FDR.component = mclapply(1:100, function(r){
  
  one.res = do.call(rbind.data.frame, lapply(method.dir, function(m){
    one.method = readRDS(paste0(dir, "/", m$dir, "/res", "_n", r, ".rds"))
    one.res = do.call(rbind.data.frame, list(
      do.call(rbind.data.frame, lapply(seq(1:length(ind.list))[!grepl("offset", names(ind.list))], function(a){
        data.frame(r = r, method = m$name, FDR.type = "FDR.offset", 
                   type = gsub("ind\\.()", "\\1", names(ind.list)[a]), measure = sum(one.method$ind.Mesor[ind.list[[a]]])/sum(one.method$ind.Mesor))
      })), 
      do.call(rbind.data.frame, lapply(seq(1:length(ind.list))[!grepl("A", names(ind.list))], function(a){
        data.frame(r = r, method = m$name, FDR.type = "FDR.A", 
                   type = gsub("ind\\.()", "\\1", names(ind.list)[a]), measure = sum(one.method$ind.amp[ind.list[[a]]])/sum(one.method$ind.amp))
      })), 
      do.call(rbind.data.frame, lapply(seq(1:length(ind.list))[!grepl("phase", names(ind.list))], function(a){
        data.frame(r = r, method = m$name, FDR.type = "FDR.phase", 
                   type = gsub("ind\\.()", "\\1", names(ind.list)[a]), measure = sum(one.method$ind.phase[ind.list[[a]]])/sum(one.method$ind.phase))
      }))
    ))

  }))
}, mc.cores = 20)
FDR.component.tab = do.call(rbind.data.frame, FDR.component) 

p = ggplot(data = FDR.component.tab, aes(x = type, y = measure, fill = method))+
  geom_boxplot()+
  facet_grid(.~FDR.type)
pdf("FDRcomponent.pdf", width = 25, height = 10)
print(p)
dev.off()