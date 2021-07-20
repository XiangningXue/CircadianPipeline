##differential rhythmicity: global test
#compare the circacompare with permutation usnig simulated data
#speed3
# Set directory and load libraries-----------------------------------------------------------
rm(list = ls())
dir = "/home/xix66/circadian/ThePipeline/diffRhyth_global_bonferroni_noR2change"
setwd(dir)

#the compare code
res.circacompare = readRDS(paste0(dir, "/res_circacompare/res", "_n", 1, ".rds"))
rthmicChange = readRDS(paste0(dir, "/res_cosinor_permutation/res", "_n", 1, ".rds"))


pdf(paste0(dir, "/Compare_Permutation_Circacompare.pdf"), width = 12, height = 36)
par(mfrow = c(9, 3))
plot(x = -log10(rthmicChange$intshiftPvalue[1:650]), y = -log10(res.circacompare$p.Mesor[1:650]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset change in 'none' genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[1:650]), y = -log10(res.circacompare$p.amp[1:650]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp change in 'none' genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[1:650]), y = -log10(res.circacompare$p.amp[1:650]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase change in 'none' genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[651:700]), y = -log10(res.circacompare$p.Mesor[651:700]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[651:700]), y = -log10(res.circacompare$p.amp[651:700]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[651:700]), y = -log10(res.circacompare$p.amp[651:700]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'offset' changed genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[701:750]), y = -log10(res.circacompare$p.Mesor[701:750]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'A' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[701:750]), y = -log10(res.circacompare$p.amp[701:750]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'A' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[701:750]), y = -log10(res.circacompare$p.amp[701:750]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'A' changed genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[751:800]), y = -log10(res.circacompare$p.Mesor[751:800]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'phase' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[751:800]), y = -log10(res.circacompare$p.amp[751:800]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'phase' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[751:800]), y = -log10(res.circacompare$p.amp[751:800]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'phase' changed genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[801:850]), y = -log10(res.circacompare$p.Mesor[801:850]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'A.phase' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[801:850]), y = -log10(res.circacompare$p.amp[801:850]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'A.phase' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[801:850]), y = -log10(res.circacompare$p.amp[801:850]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'A.phase' changed genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[851:900]), y = -log10(res.circacompare$p.Mesor[851:900]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'A.offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[851:900]), y = -log10(res.circacompare$p.amp[851:900]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'A.offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[851:900]), y = -log10(res.circacompare$p.amp[851:900]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'A.offset' changed genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[901:950]), y = -log10(res.circacompare$p.Mesor[901:950]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'phase.offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[901:950]), y = -log10(res.circacompare$p.amp[901:950]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'phase.offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[901:950]), y = -log10(res.circacompare$p.amp[901:950]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'phase.offset' changed genes")
curve((x), col = "green", add = TRUE)

plot(x = -log10(rthmicChange$intshiftPvalue[951:1000]), y = -log10(res.circacompare$p.Mesor[951:1000]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of offset  in 'A.phase.offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$AshiftPvalue[951:1000]), y = -log10(res.circacompare$p.amp[951:1000]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of Amp in 'A.phase.offset' changed genes")
curve((x), col = "green", add = TRUE)
plot(x = -log10(rthmicChange$phasePvalue[951:1000]), y = -log10(res.circacompare$p.amp[951:1000]), xlab = "permutation", ylab = "circacompare", main = "-log10Pval of phase in 'A.phase.offset' changed genes")
curve((x), col = "green", add = TRUE)

dev.off()