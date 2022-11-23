library(ape)
library(picante)
library(tidyverse)
library(phytools)
library(caper)
library(ggrepel)
library(cowplot)
library(ggsci)
library(patchwork)

modPgls.SEy = function (model, data, corClass = corBrownian, tree, se = NULL,
                        method = c("REML", "ML"), interval = c(0, 1000), corClassValue=1, sig2e=NULL, ...)
{
  Call <- match.call()
  corfunc <- corClass
  spp <- rownames(data)
  data <- cbind(data, spp)
  if (is.null(se))
    se <- setNames(rep(0, Ntip(tree)), tree$tip.label)[spp]
  else se <- se[spp]

  lk <- function(sig2e, data, tree, model, ve, corfunc, spp) {
    tree$edge.length <- tree$edge.length * sig2e
    ii <- sapply(1:Ntip(tree), function(x, e) which(e ==
                                                      x), e = tree$edge[, 2])
    tree$edge.length[ii] <- tree$edge.length[ii] + ve[tree$tip.label]
    vf <- diag(vcv(tree))[spp]
    w <- varFixed(~vf)
    COR <- corfunc(corClassValue, tree, form = ~spp, ...)
    fit <- gls(model, data = cbind(data, vf), correlation = COR,
               method = method, weights = w)
    -logLik(fit)
  }

  if (is.null(sig2e)) {
    fit <- optimize(lk, interval = interval, data = data, tree = tree,
                    model = model, ve = se^2, corfunc = corfunc, spp = spp)
    sig2e=fit$minimum
  }

  tree$edge.length <- tree$edge.length * sig2e
  ii <- sapply(1:Ntip(tree), function(x, e) which(e == x),
               e = tree$edge[, 2])
  tree$edge.length[ii] <- tree$edge.length[ii] + se[tree$tip.label]^2
  vf <- diag(vcv(tree))[spp]
  w <- varFixed(~vf)
  obj <- gls(model, data = cbind(data, vf), correlation = corfunc(corClassValue,
                                                                  tree, form = ~spp, ...), weights = w, method = method)
  obj$call <- Call
  obj$sig2e <- sig2e
  obj
}

#Internal function
pglsSEyPagelToOptimizeLambda=function(lambda,model,data,tree,...) {
  -logLik(modPgls.SEy(model=model,data=data,tree=tree,corClassValue=lambda,corClass=corPagel,fixed=T,...)) #Returns -logLikelihood of the pgls.SEy model with lambda fixed to the value of the lambda argument. sig2e will be optimized within modPgls.SEy unless given as an argument here
}

#Function intended for users
pglsSEyPagel=function(model, data, tree, lambdaInterval=c(0,1),...){
  optimizedModel=optimize(pglsSEyPagelToOptimizeLambda,lambdaInterval,model=model,data=data,tree=tree,...) #Optimizes lambda in the lambdaInterval using the pglsSEyPagelToOptimizeLambda function
  return(modPgls.SEy(model=model,data=data,tree=tree,corClass=corPagel,fixed=T,corClassValue=optimizedModel$minimum,...)) #Returns the final model fit
}

## Here I split the csv into class but do whatever you want
Data<-read.csv(file="min20RAD_UPDATE.csv")
View(Data)

tree<-read.tree(file="min10radtree.nwk")



Data$Species <- gsub(" ", "_", Data$Species)
includedSpecies <- Data$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
length(pruned.tree$tip.label)
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$Species %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]
rownames(Data)<-Data$Species
SE<-setNames(Data$SE,Data$Species)[rownames(Data)]
Data<-na.omit(Data)
##Model
AUC.10.GY <- pglsSEyPagel(NeoplasiaPrevalence~log10(AUC10Gy), data=Data,
                            tree=pruned.tree,method="ML",se=SE)
summary(AUC.10.GY)
auc<-lm(NeoplasiaPrevalence~log10(AUC10Gy), data=Data)
summary(auc)$r.squared

r.AUC.10.GY <- summary(AUC.10.GY)$corBeta
r.AUC.10.GY <- format(r.AUC.10.GY[2,1])
r.AUC.10.GY <-signif(as.numeric(r.AUC.10.GY)^2, digits = 2)
ld.v.AUC.10.GY<- summary(AUC.10.GY)$modelStruct$corStruct
ld.v.AUC.10.GY <- signif(ld.v.AUC.10.GY[1], digits =2)
p.v.AUC.10.GY<-summary(AUC.10.GY)$tTable
p.v.AUC.10.GY<-signif(p.v.AUC.10.GY[2,4], digits = 2)



ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(AUC10Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(AUC.10.GY)[1]*100, slope =  coef(AUC.10.GY)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) Cell Count AUC % of Untreated Control \n 10Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC 10Gy Radiation",
        subtitle =bquote(p-value:.(p.v.AUC.10.GY)~R^2:.(r.AUC.10.GY)~Lambda:.(ld.v.AUC.10.GY)))+
  guides(col=FALSE)

ggsave(filename='RADAUC10.png', width=9.5, height=7, limitsize=FALSE,bg="white")

#AUC 2



##Model
AUC.2.GY <- pglsSEyPagel(NeoplasiaPrevalence~log10(AUC2Gy), data=Data,
                          tree=pruned.tree,method="ML",se=SE)

r.AUC.2.GY <- summary(AUC.2.GY)$corBeta
r.AUC.2.GY <- format(r.AUC.2.GY[2,1])
r.AUC.2.GY <-signif(as.numeric(r.AUC.2.GY)^2, digits = 2)
ld.v.AUC.2.GY<- summary(AUC.2.GY)$modelStruct$corStruct
ld.v.AUC.2.GY <- signif(ld.v.AUC.2.GY[1], digits = 2)
p.v.AUC.2.GY<-summary(AUC.2.GY)$tTable
p.v.AUC.2.GY<-signif(p.v.AUC.2.GY[2,4], digits = 2)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(AUC2Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(AUC.2.GY)[1]*100, slope =  coef(AUC.2.GY)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) Cell Count AUC % of Untreated Control \n 2Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC 2GY Radiation",
        subtitle =bquote(p-value:.(p.v.AUC.2.GY)~R^2:.(r.AUC.2.GY)~Lambda:.(ld.v.AUC.2.GY)))+
  guides(col=FALSE)

ggsave(filename='RADAUC2.png', width=9.5, height=7, limitsize=FALSE,bg="white")

#AUC .4

##Model
AUC.04.GY <- pglsSEyPagel(NeoplasiaPrevalence~log10(AUC04Gy), data=Data,
                         tree=pruned.tree,method="ML",se=SE)

r.AUC.04.GY <- summary(AUC.04.GY)$corBeta
r.AUC.04.GY <- format(r.AUC.04.GY[2,1])
r.AUC.04.GY <-signif(as.numeric(r.AUC.04.GY)^2, digits = 2)
ld.v.AUC.04.GY<- summary(AUC.04.GY)$modelStruct$corStruct
ld.v.AUC.04.GY <- signif(ld.v.AUC.04.GY[1], digits = 2)
p.v.AUC.04.GY<-summary(AUC.04.GY)$tTable
p.v.AUC.04.GY<-signif(p.v.AUC.04.GY[2,4], digits = 2)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(AUC04Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(AUC.04.GY)[1]*100, slope =  coef(AUC.04.GY)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) Cell Count AUC of Untreated Control \n 0.4Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC 0.4Gy Radiation",
        subtitle =bquote(p-value:.(p.v.AUC.04.GY)~R^2:.(r.AUC.04.GY)~Lambda:.(ld.v.AUC.04.GY)))+
  guides(col=FALSE)

ggsave(filename='RADAUC04.png', width=9.5, height=7, limitsize=FALSE,bg="white")

##CHANGE 


#10GY Change
Change.10.GY <- pglsSEyPagel(NeoplasiaPrevalence~log10(Change10Gy), data=Data,
                          tree=pruned.tree,method="ML",se=SE)

r.Change.10.GY <- summary(Change.10.GY)$corBeta
r.Change.10.GY <- format(r.Change.10.GY[2,1])
r.Change.10.GY <-signif(as.numeric(r.Change.10.GY)^2, digits = 2)
ld.v.Change.10.GY<- summary(Change.10.GY)$modelStruct$corStruct
ld.v.Change.10.GY <- signif(ld.v.Change.10.GY[1], digits =2)
p.v.Change.10.GY<-summary(Change.10.GY)$tTable
p.v.Change.10.GY<-signif(p.v.Change.10.GY[2,4], digits = 2)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(Change10Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(Change.10.GY)[1]*100, slope =  coef(Change.10.GY)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("AUC Cell Count % Area Change from UT 10 Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC Area Change From UT 10Gy Radiation",
       subtitle =bquote(p-value:.(p.v.Change.10.GY)~R^2:.(r.Change.10.GY)~Lambda:.(ld.v.Change.10.GY)))+
  guides(col=FALSE)

ggsave(filename='RADChange10.png', width=9.5, height=7, limitsize=FALSE,bg="white")

#Change 2



##Model
Change.2.GY <- pglsSEyPagel(NeoplasiaPrevalence~log10(Change2Gy), data=Data,
                         tree=pruned.tree,method="ML",se=SE)

r.Change.2.GY <- summary(Change.2.GY)$corBeta
r.Change.2.GY <- format(r.Change.2.GY[2,1])
r.Change.2.GY <-signif(as.numeric(r.Change.2.GY)^2, digits = 2)
ld.v.Change.2.GY<- summary(Change.2.GY)$modelStruct$corStruct
ld.v.Change.2.GY <- signif(ld.v.Change.2.GY[1], digits = 2)
p.v.Change.2.GY<-summary(Change.2.GY)$tTable
p.v.Change.2.GY<-signif(p.v.Change.2.GY[2,4], digits = 2)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(Change2Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(Change.2.GY)[1]*100, slope =  coef(Change.2.GY)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("AUC Cell Count % Area Change from UT 2 Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC Area Change From UT 2Gy Radiation",
       subtitle =bquote(p-value:.(p.v.Change.2.GY)~R^2:.(r.Change.2.GY)~Lambda:.(ld.v.Change.2.GY)))+
  guides(col=FALSE)

ggsave(filename='RADChange2.png', width=9.5, height=7, limitsize=FALSE,bg="white")

#Change .4

##Model
Change.04.GY <- pglsSEyPagel(NeoplasiaPrevalence~log10(Change04Gy), data=Data,
                          tree=pruned.tree,method="ML",se=SE)

r.Change.04.GY <- summary(Change.04.GY)$corBeta
r.Change.04.GY <- format(r.Change.04.GY[2,1])
r.Change.04.GY <-signif(as.numeric(r.Change.04.GY)^2, digits = 3)
ld.v.Change.04.GY<- summary(Change.04.GY)$modelStruct$corStruct
ld.v.Change.04.GY <- signif(ld.v.Change.04.GY[1], digits = 2)
p.v.Change.04.GY<-summary(Change.04.GY)$tTable
p.v.Change.04.GY<-signif(p.v.Change.04.GY[2,4], digits = 2)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(Change04Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(Change.04.GY)[1]*100, slope =  coef(Change.04.GY)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("AUC Cell Count % Area Change from UT .4 Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC Area Change From UT .4Gy Radiation",
       subtitle =bquote(p-value:.(p.v.Change.04.GY)~R^2:.(r.Change.04.GY)~Lambda:.(ld.v.Change.04.GY)))+
  guides(col=FALSE)

ggsave(filename='RADChange04.png', width=9.5, height=7, limitsize=FALSE,bg="white")





#malignancy


#10GY Change
Change.10.GY.mal <- pglsSEyPagel(MalignancyPrevalence~log10(Change10Gy), data=Data,
                             tree=pruned.tree,method="ML",se=SE)

r.Change.10.GY.mal <- summary(Change.10.GY.mal)$corBeta
r.Change.10.GY.mal <- format(r.Change.10.GY.mal[2,1])
r.Change.10.GY.mal <-signif(as.numeric(r.Change.10.GY.mal)^2, digits = 2)
ld.v.Change.10.GY.mal<- summary(Change.10.GY.mal)$modelStruct$corStruct
ld.v.Change.10.GY.mal <- signif(ld.v.Change.10.GY.mal[1], digits =2)
p.v.Change.10.GY.mal<-summary(Change.10.GY.mal)$tTable
p.v.Change.10.GY.mal<-signif(p.v.Change.10.GY.mal[2,4], digits = 2)


ggplot(Data, aes(y=MalignancyPrevalence*100, x=log10(Change10Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(Change.10.GY.mal)[1]*100, slope =  coef(Change.10.GY.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("AUC Cell Count % Area Change from UT 10 Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Malignancy vs. AUC Area Change From UT 10Gy Radiation",
       subtitle =bquote(p-value:.(p.v.Change.10.GY.mal)~R^2:.(r.Change.10.GY.mal)~Lambda:.(ld.v.Change.10.GY.mal)))+
  guides(col=FALSE)

ggsave(filename='RADChange10mal.png', width=9.5, height=7, limitsize=FALSE,bg="white")

#Change 2



##Model
Change.2.GY.mal <- pglsSEyPagel(MalignancyPrevalence~log10(Change2Gy), data=Data,
                            tree=pruned.tree,method="ML",se=SE)

r.Change.2.GY.mal <- summary(Change.2.GY.mal)$corBeta
r.Change.2.GY.mal <- format(r.Change.2.GY.mal[2,1])
r.Change.2.GY.mal <-signif(as.numeric(r.Change.2.GY.mal)^2, digits = 2)
ld.v.Change.2.GY.mal<- summary(Change.2.GY.mal)$modelStruct$corStruct
ld.v.Change.2.GY.mal <- signif(ld.v.Change.2.GY.mal[1], digits = 2)
p.v.Change.2.GY.mal<-summary(Change.2.GY.mal)$tTable
p.v.Change.2.GY.mal<-signif(p.v.Change.2.GY.mal[2,4], digits = 2)


ggplot(Data, aes(y=MalignancyPrevalence*100, x=log10(Change2Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(Change.2.GY.mal)[1]*100, slope =  coef(Change.2.GY.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("AUC Cell Count % Area Change from UT 2 Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Malignancy vs. AUC Area Change From UT 2Gy Radiation",
       subtitle =bquote(p-value:.(p.v.Change.2.GY.mal)~R^2:.(r.Change.2.GY.mal)~Lambda:.(ld.v.Change.2.GY.mal)))+
  guides(col=FALSE)

ggsave(filename='RADChange2mal.png', width=9.5, height=7, limitsize=FALSE,bg="white")

#Change .4

##Model
Change.04.GY.mal <- pglsSEyPagel(MalignancyPrevalence~log10(Change04Gy), data=Data,
                             tree=pruned.tree,method="ML",se=SE)

r.Change.04.GY.mal <- summary(Change.04.GY.mal)$corBeta
r.Change.04.GY.mal <- format(r.Change.04.GY.mal[2,1])
r.Change.04.GY.mal <-signif(as.numeric(r.Change.04.GY.mal)^2, digits = 2)
ld.v.Change.04.GY.mal<- summary(Change.04.GY.mal)$modelStruct$corStruct
ld.v.Change.04.GY.mal <- signif(ld.v.Change.04.GY.mal[1], digits = 2)
p.v.Change.04.GY.mal<-summary(Change.04.GY.mal)$tTable
p.v.Change.04.GY.mal<-signif(p.v.Change.04.GY.mal[2,4], digits = 2)


ggplot(Data, aes(y=MalignancyPrevalence*100, x=log10(Change04Gy))) +
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(Change.04.GY.mal)[1]*100, slope =  coef(Change.04.GY.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("AUC Cell Count % Area Change from UT .4 Gy Radiation") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=common_name))+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,150,250,394),
             labels =  c(20,150,250,394))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Malignancy vs. AUC Area Change From UT .4Gy Radiation",
       subtitle =bquote(p-value:.(p.v.Change.04.GY.mal)~R^2:.(r.Change.04.GY.mal)~Lambda:.(ld.v.Change.04.GY.mal)))+
  guides(col=FALSE)

ggsave(filename='RADChange04mal.png', width=9.5, height=7, limitsize=FALSE,bg="white")
