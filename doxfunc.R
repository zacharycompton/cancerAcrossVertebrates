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
Data<-read.csv(file="min20DOX.csv")
View(Data)

tree<-read.tree(file="min20Fixed516.nwk")
length(tree$tip.label)

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
##Model
ANVFold72.1 <- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVAUCFold72.1), data=Data,
                            tree=pruned.tree,method="ML",se=SE)

r.v.ANVFold72.1 <- summary(ANVFold72.1)$corBeta
r.v.ANVFold72.1 <- format(r.v.ANVFold72.1[2,1])
r.v.ANVFold72.1 <-signif(as.numeric(r.v.ANVFold72.1)^2, digits= 3)
ld.v.ANVFold72.1<- summary(ANVFold72.1)$modelStruct$corStruct
ld.v.ANVFold72.1 <- signif(ld.v.ANVFold72.1[1])
p.v.ANVFold72.1<-summary(ANVFold72.1)$tTable
p.v.ANVFold72.1<-signif(p.v.ANVFold72.1[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(ANVAUCFold72.1))) +
  scale_x_continuous(
    limits = c(-.01,1.397),
    breaks = c(-.01,1,1.397),
    labels = c(0.97,10,25)
  )+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVFold72.1)[1]*100, slope =  coef(ANVFold72.1)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) AUC AnV Fold \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=Data$common_name, size = 300))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(title = "Neoplasia vs. AUC 1uM Doxorubicin",
       subtitle =bquote(p-value:.(p.v.ANVFold72.1)~R^2:.(r.v.ANVFold72.1)~Lambda:.(ld.v.ANVFold72.1)))+

  guides(col=FALSE)

ggsave(filename='anvfold721.png', width=13, height=10, limitsize=FALSE,bg="white")


#fold 72.33
##Model
ANVFold72.33<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVAUCFold72.33), data=Data,
                            tree=pruned.tree,method="ML",se=SE)

r.v.ANVFold72.33 <- summary(ANVFold72.33)$corBeta
r.v.ANVFold72.33 <- format(r.v.ANVFold72.33[2,1])
r.v.ANVFold72.33 <-signif(as.numeric(r.v.ANVFold72.33)^2, digits= 3)
ld.v.ANVFold72.33<- summary(ANVFold72.33)$modelStruct$corStruct
ld.v.ANVFold72.33<- signif(ld.v.ANVFold72.33[1])
p.v.ANVFold72.33<-summary(ANVFold72.33)$tTable
p.v.ANVFold72.33<-signif(p.v.ANVFold72.33[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(ANVAUCFold72.33))) +
  scale_x_continuous(
    limits = c(-.3,1),
    breaks = c(-.3,.69897,1),
    labels = c(0.5,5,10)
  )+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVFold72.33)[1]*100, slope =  coef(ANVFold72.33)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) AUC AnV Fold \n % Increase to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
   guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom") +
  labs(title = "Neoplasia vs. AUC 0.33uM Doxorubicin",
       subtitle =bquote(p-value:.(p.v.ANVFold72.33)~R^2:.(r.v.ANVFold72.33)~Lambda:.(ld.v.ANVFold72.33)))+
  guides(col=FALSE)

ggsave(filename='anvfold72.33.png', width=13, height=10, limitsize=FALSE,bg="white")


##Model
ANVFold72.11<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVAUCFold72.11), data=Data,
                            tree=pruned.tree,method="ML",se=SE)

r.v.ANVFold72.11 <- summary(ANVFold72.11)$corBeta
r.v.ANVFold72.11<- format(r.v.ANVFold72.11[2,1])
r.v.ANVFold72.11 <-signif(as.numeric(r.v.ANVFold72.11)^2, digits= 3)
ld.v.ANVFold72.11<- summary(ANVFold72.11)$modelStruct$corStruct
ld.v.ANVFold72.11<- signif(ld.v.ANVFold72.11[1], digits = 3)
p.v.ANVFold72.11<-summary(ANVFold72.11)$tTable
p.v.ANVFold72.11<-signif(p.v.ANVFold72.11[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(ANVAUCFold72.11))) +
  scale_x_continuous(
    limits = c(-.02,.69897),
    breaks = c(-.02,.69897),
    labels = c(0.95,5)
  )+
  scale_color_manual(values=c("#631879FF"))+

  geom_abline(intercept = coef(ANVFold72.11)[1]*100, slope =  coef(ANVFold72.11)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) AUC AnV Fold \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom") +
  labs(title = "Neoplasia vs. AUC 0.11uM Doxorubicin",
       subtitle =bquote(p-value:.(p.v.ANVFold72.11)~R^2:.(r.v.ANVFold72.11)~Lambda:.(ld.v.ANVFold72.11)))+
  guides(col=FALSE)

ggsave(filename='anvfold72.11.png', width=13, height=10, limitsize=FALSE,bg="white")



#cell death 72.1
##Model
ANVDeath72.1 <- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVCellDeath72.1), data=Data,
                             tree=pruned.tree,method="ML",se=SE)

r.v.ANVDeath72.1 <- summary(ANVDeath72.1)$corBeta
r.v.ANVDeath72.1 <- format(r.v.ANVDeath72.1[2,1])
r.v.ANVDeath72.1 <-signif(as.numeric(r.v.ANVDeath72.1)^2, digits= 3)
ld.v.ANVDeath72.1<- summary(ANVDeath72.1)$modelStruct$corStruct
ld.v.ANVDeath72.1 <- signif(ld.v.ANVDeath72.1[1])
p.v.ANVDeath72.1<-summary(ANVDeath72.1)$tTable
p.v.ANVDeath72.1<-signif(p.v.ANVDeath72.1[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(ANVCellDeath72.1))) +
  scale_x_continuous(
    limits = c(-.4,1.74),
    breaks = c(-.4,1,1.47712,1.74),
    labels = c(0.4,10,30,60)
  )+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVDeath72.1)[1]*100, slope =  coef(ANVDeath72.1)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) % Cell Death at 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom") +
  labs(title = "Neoplasia vs. Cell Death 1uM Doxorubicin",
       subtitle =bquote(p-value:.(p.v.ANVDeath72.1)~R^2:.(r.v.ANVDeath72.1)~Lambda:.(ANVDeath72.1)))+
  guides(col=FALSE)

ggsave(filename='celldeath72.1.png', width=13, height=10, limitsize=FALSE,bg="white")


#cell death 72.33


##Model
ANVDeath72.33<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVCellDeath72.33), data=Data,
                             tree=pruned.tree,method="ML",se=SE)

r.v.ANVDeath72.33 <- summary(ANVDeath72.33)$corBeta
r.v.ANVDeath72.33 <- format(r.v.ANVDeath72.33[2,1])
r.v.ANVDeath72.33 <-signif(as.numeric(r.v.ANVDeath72.33)^2, digits= 3)
ld.v.ANVDeath72.33<- summary(ANVDeath72.33)$modelStruct$corStruct
ld.v.ANVDeath72.33<- signif(ld.v.ANVDeath72.33[1])
p.v.ANVDeath72.33<-summary(ANVDeath72.33)$tTable
p.v.ANVDeath72.33<-signif(p.v.ANVDeath72.33[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(ANVCellDeath72.33))) +
  scale_x_continuous(
    limits = c(-.8,1.74),
    breaks = c(-.8,1,1.47712,1.74),
    labels = c(.16,10,30,60)
  )+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVDeath72.33)[1]*100, slope =  coef(ANVDeath72.33)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("% Cell Death at 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom") +
  labs(title = "Neoplasia vs. Cell Death 0.33uM Doxorubicin",
       subtitle =bquote(p-value:.(p.v.ANVDeath72.33)~R^2:.(r.v.ANVDeath72.33)~Lambda:.(ANVDeath72.33)))+
  guides(col=FALSE)

ggsave(filename='celldeath72.33.png', width=13, height=10, limitsize=FALSE,bg="white")

#cell death 72.11

##Model
ANVDeath72.11<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVCellDeath72.11), data=Data,
                             tree=pruned.tree,method="ML",se=SE)

r.v.ANVDeath72.11 <- summary(ANVDeath72.11)$corBeta
r.v.ANVDeath72.11<- format(r.v.ANVDeath72.11[2,1])
r.v.ANVDeath72.11 <-signif(as.numeric(r.v.ANVDeath72.11)^2, digits= 3)
ld.v.ANVDeath72.11<- summary(ANVDeath72.11)$modelStruct$corStruct
ld.v.ANVDeath72.11<- signif(ld.v.ANVDeath72.11[1], digits = 3)
p.v.ANVDeath72.11<-summary(ANVDeath72.11)$tTable
p.v.ANVDeath72.11<-signif(p.v.ANVDeath72.11[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=log10(ANVCellDeath72.11))) +
  scale_x_continuous(
    limits = c(-.8,1.47712),
    breaks = c(-.8,1,1.47712),
    labels = c(.16,10,30)
  )+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVDeath72.11)[1]*100, slope =  coef(ANVDeath72.11)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) % Cell Death at 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom") +
  labs(title = "Neoplasia vs. Cell Death 0.11uM Doxorubicin",
       subtitle =bquote(p-value:.(p.v.ANVDeath72.11)~R^2:.(r.v.ANVDeath72.11)~Lambda:.(ANVDeath72.11)))+
  guides(col=FALSE)

ggsave(filename='celldeath72.11.png', width=13, height=10, limitsize=FALSE,bg="white")


