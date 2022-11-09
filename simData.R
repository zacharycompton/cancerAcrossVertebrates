tree<-pbtree(n=97,scale=1)
## simulate data under model="BM"
x<-fastBM(tree)
plotTree(tree)
x<-fastBM(rescale(tree,model="OU",0.4))




fitEB_sim <- fitContinuous(tree, x, model="EB")
fitBM_sim <- fitContinuous(tree, x, model="BM")
fitOU_sim <- fitContinuous(tree, x, model="OU")

aic_sim <- setNames(c(AIC(fitOUMCMC_sim), AIC(fitBMMCMC_sim),
                      AIC(fitEBMCMC_sim)), c("OU", "BM", "EB"))
aic_sim

fitOUMCMC_sim <- fitContinuousMCMC(tree, x, model="SSP", Ngens = 10000,
                                   sampleFreq = 1000, printFreq = 1000)
fitBMMCMC_sim <- fitContinuousMCMC(tree, x, model="BM", Ngens = 10000,
                                   sampleFreq = 1000, printFreq = 1000)
fitEBMCMC_sim <- fitContinuousMCMC(tree, x, model="ACDC.exp", Ngens = 10000,
                                   sampleFreq = 1000, printFreq = 1000)
