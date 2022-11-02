tree<-pbtree(n=98,scale=1)
## simulate data under model="BM"
x<-fastBM(tree)

x<-fastBM(rescale(tree,model="EB",0.4))

library(geiger)



fitEB_sim <- fitContinuous(tree, x, model="EB")
fitBM_sim <- fitContinuous(tree, x, model="BM")
fitOU_sim <- fitContinuous(tree, x, model="OU")

aic_sim <- setNames(c(AIC(fitBM_sim), AIC(fitEB_sim),
                      AIC(fitOU_sim)), c("BM", "EB", "OU"))
aic_sim
