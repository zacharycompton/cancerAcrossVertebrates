library(ape)
library(picante)
library(tidyverse)
library(phytools)
library(caper)
library(ggrepel)
library(cowplot)
library(ggsci)

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
Data <- read.csv("min50Functional.csv")
Data <- Data[-c(9), ]
View(Data)
Data <- mutate(Data, SE = sqrt(1/(Data$TotalRecords)))
Data <- Data[!(Data$TotalRecords<20),]
Data[Data==-1]<-NA
Data[Data < 0] <-NA
#Subset if you only want a certain Species
#Data <- Data[ which(Species == "Mammalia"),]

tree <- read.tree("min20Fixed516.nwk")


#prune the tree to match the data
Data$Species <- gsub(" ", "_", Data$Species)
includedSpecies <- Data$Species
#pruning the tree
tree$tip.label <- newtips
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
#Removing discrepencies
Data$Keep <- Data$Species %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]

##Clean up for PGLS.sey 
cutData <- Data[,c(2,3,4,5,6,9),drop=FALSE] 
View(cutData)
## species labels as row names
rownames(cutData)<-cutData$Species
## pull out the SEs
SE<-setNames(cutData$SE,cutData$Species)[rownames(cutData)]


##Model
ANVFold72.1 <- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVAUCFold72.1), data=cutData,
                            tree=pruned.tree,method="ML",se=SE)

r.v.ANVFold72.1 <- summary(ANVFold72.1)$corBeta
r.v.ANVFold72.1 <- format(r.v.ANVFold72.1[2,1])
r.v.ANVFold72.1 <-signif(as.numeric(r.v.ANVFold72.1)^2, digits= 3)
ld.v.ANVFold72.1<- summary(ANVFold72.1)$modelStruct$corStruct
ld.v.ANVFold72.1 <- signif(ld.v.ANVFold72.1[1])
p.v.ANVFold72.1<-summary(ANVFold72.1)$tTable
p.v.ANVFold72.1<-signif(p.v.ANVFold72.1[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=(ANVAUCFold72.1))) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVFold72.1)[1]*100, slope =  coef(ANVFold72.1)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log 10) AUC AnV Fold \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=cutData$common_name, size = 300))+
  scale_size(name   = "Total Necropsies",
             breaks = c(55,150,250,350),
             labels =  c(55,150,250,350))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+ 
  guides(col=FALSE)

ggsave(filename='anvfold721.png', width=13, height=10, limitsize=FALSE,bg="white")


#fold 72.33
cutData <- Data[,c(2,3,4,5,6,10),drop=FALSE] 
View(cutData)
## species labels as row names
rownames(cutData)<-cutData$Species
## pull out the SEs
SE<-setNames(cutData$SE,cutData$Species)[rownames(cutData)]


##Model
ANVFold72.33<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVAUCFold72.33), data=cutData,
                            tree=pruned.tree,method="ML",se=SE)

r.v.ANVFold72.33 <- summary(ANVFold72.33)$corBeta
r.v.ANVFold72.33 <- format(r.v.ANVFold72.33[2,1])
r.v.ANVFold72.33 <-signif(as.numeric(r.v.ANVFold72.33)^2, digits= 3)
ld.v.ANVFold72.33<- summary(ANVFold72.33)$modelStruct$corStruct
ld.v.ANVFold72.33<- signif(ld.v.ANVFold72.33[1])
p.v.ANVFold72.33<-summary(ANVFold72.33)$tTable
p.v.ANVFold72.33<-signif(p.v.ANVFold72.33[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=(ANVAUCFold72.33))) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVFold72.33)[1]*100, slope =  coef(ANVFold72.33)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Log10 AUC AnV Fold \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  labs(title = "Neoplasia Prevalence vs. AUC AnV Fold \nIncrease to 72hr .33uM in Mammals",  
       subtitle = (paste("p-value:", p.v.ANVFold72.33,"  ","R^2:",r.v.ANVFold72.33,"  ","Lambda:",ld.v.ANVFold72.33))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")

#fold 72.11
cutData <- Data[,c(2,3,4,5,6,11),drop=FALSE] 
View(cutData)
## species labels as row names
rownames(cutData)<-cutData$Species
## pull out the SEs
SE<-setNames(cutData$SE,cutData$Species)[rownames(cutData)]


##Model
ANVFold72.11<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVAUCFold72.11), data=cutData,
                            tree=pruned.tree,method="ML",se=SE)

r.v.ANVFold72.11 <- summary(ANVFold72.11)$corBeta
r.v.ANVFold72.11<- format(r.v.ANVFold72.11[2,1])
r.v.ANVFold72.11 <-signif(as.numeric(r.v.ANVFold72.11)^2, digits= 3)
ld.v.ANVFold72.11<- summary(ANVFold72.11)$modelStruct$corStruct
ld.v.ANVFold72.11<- signif(ld.v.ANVFold72.11[1], digits = 3)
p.v.ANVFold72.11<-summary(ANVFold72.11)$tTable
p.v.ANVFold72.11<-signif(p.v.ANVFold72.11[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=(ANVAUCFold72.11))) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual(values=c("#631879FF"))+
  
  geom_abline(intercept = coef(ANVFold72.11)[1]*100, slope =  coef(ANVFold72.11)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Log10 AUC AnV Fold \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  #geom_smooth(method=lm,colour = "grey" ,se=FALSE, size=1.5, ) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  labs(title = "Neoplasia Prevalence vs. AUC AnV Fold \nIncrease to 72hr .11uM in Mammals",  
       subtitle = (paste("p-value:", p.v.ANVFold72.11,"  ","R^2:",r.v.ANVFold72.11,"  ","Lambda:",ld.v.ANVFold72.11))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")





#cell death 72.1

cutData <- Data[,c(2,3,4,5,6,12),drop=FALSE] 
View(cutData)
## species labels as row names
rownames(cutData)<-cutData$Species
## pull out the SEs
SE<-setNames(cutData$SE,cutData$Species)[rownames(cutData)]


##Model
ANVDeath72.1 <- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVCellDeath72.1), data=cutData,
                             tree=pruned.tree,method="ML",se=SE)

r.v.ANVDeath72.1 <- summary(ANVDeath72.1)$corBeta
r.v.ANVDeath72.1 <- format(r.v.ANVDeath72.1[2,1])
r.v.ANVDeath72.1 <-signif(as.numeric(r.v.ANVDeath72.1)^2, digits= 3)
ld.v.ANVDeath72.1<- summary(ANVDeath72.1)$modelStruct$corStruct
ld.v.ANVDeath72.1 <- signif(ld.v.ANVDeath72.1[1])
p.v.ANVDeath72.1<-summary(ANVDeath72.1)$tTable
p.v.ANVDeath72.1<-signif(p.v.ANVDeath72.1[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=(ANVCellDeath72.1))) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVDeath72.1)[1]*100, slope =  coef(ANVDeath72.1)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Log10 AUC AnV %Cell Death \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  #geom_smooth(method=lm,colour = "grey" ,se=FALSE, size=1.5, ) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  labs(title = "Neoplasia Prevalence vs. AUC AnV Cell Death \nIncrease to 72hr 1uM in Mammals",  
       subtitle = (paste("p-value:", p.v.ANVDeath72.1,"  ","R^2:",r.v.ANVDeath72.1,"  ","Lambda:",ld.v.ANVDeath72.1))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")


#cell death 72.33
cutData <- Data[,c(2,3,4,5,6,13),drop=FALSE] 
View(cutData)
## species labels as row names
rownames(cutData)<-cutData$Species
## pull out the SEs
SE<-setNames(cutData$SE,cutData$Species)[rownames(cutData)]


##Model
ANVDeath72.33<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVCellDeath72.33), data=cutData,
                             tree=pruned.tree,method="ML",se=SE)

r.v.ANVDeath72.33 <- summary(ANVDeath72.33)$corBeta
r.v.ANVDeath72.33 <- format(r.v.ANVDeath72.33[2,1])
r.v.ANVDeath72.33 <-signif(as.numeric(r.v.ANVDeath72.33)^2, digits= 3)
ld.v.ANVDeath72.33<- summary(ANVDeath72.33)$modelStruct$corStruct
ld.v.ANVDeath72.33<- signif(ld.v.ANVDeath72.33[1])
p.v.ANVDeath72.33<-summary(ANVDeath72.33)$tTable
p.v.ANVDeath72.33<-signif(p.v.ANVDeath72.33[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=(ANVCellDeath72.33))) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVDeath72.33)[1]*100, slope =  coef(ANVDeath72.33)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Log10 AUC AnV %Cell Death \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  #geom_smooth(method=lm,colour = "grey" ,se=FALSE, size=1.5, ) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  labs(title = "Neoplasia Prevalence vs. AUC AnV Cell Death \nIncrease to 72hr .33uM in Mammals",  
       subtitle = (paste("p-value:", p.v.ANVDeath72.33,"  ","R^2:",r.v.ANVDeath72.33,"  ","Lambda:",ld.v.ANVDeath72.33))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")

#cell death 72.11
cutData <- Data[,c(2,3,4,5,6,14),drop=FALSE] 
View(cutData)
## species labels as row names
rownames(cutData)<-cutData$Species
## pull out the SEs
SE<-setNames(cutData$SE,cutData$Species)[rownames(cutData)]


##Model
ANVDeath72.11<- pglsSEyPagel(NeoplasiaPrevalence~log10(ANVCellDeath72.11), data=cutData,
                             tree=pruned.tree,method="ML",se=SE)

r.v.ANVDeath72.11 <- summary(ANVDeath72.11)$corBeta
r.v.ANVDeath72.11<- format(r.v.ANVDeath72.11[2,1])
r.v.ANVDeath72.11 <-signif(as.numeric(r.v.ANVDeath72.11)^2, digits= 3)
ld.v.ANVDeath72.11<- summary(ANVDeath72.11)$modelStruct$corStruct
ld.v.ANVDeath72.11<- signif(ld.v.ANVDeath72.11[1], digits = 3)
p.v.ANVDeath72.11<-summary(ANVDeath72.11)$tTable
p.v.ANVDeath72.11<-signif(p.v.ANVDeath72.11[2,4], digits = 3)


ggplot(Data, aes(y=NeoplasiaPrevalence*100, x=(ANVCellDeath72.11))) +
  scale_x_continuous(trans = 'log10')+
  scale_color_manual(values=c("#631879FF"))+
  geom_abline(intercept = coef(ANVDeath72.11)[1]*100, slope =  coef(ANVDeath72.11)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Log10 AUC AnV %Cell Death \nIncrease to 72hr") +
  geom_point(aes(colour= Keep, size = TotalRecords)) +
  #geom_smooth(method=lm,colour = "grey" ,se=FALSE, size=1.5, ) +
  geom_text_repel(aes(label=ifelse((NeoplasiaPrevalence > 0) | NeoplasiaPrevalence < 1,as.character(common_name),'')))+
  labs(title = "Neoplasia Prevalence vs. AUC AnV Cell Death \nIncrease to 72hr .11uM in Mammals",  
       subtitle = (paste("p-value:", p.v.ANVDeath72.11,"  ","R^2:",r.v.ANVDeath72.11,"  ","Lambda:",ld.v.ANVDeath72.11))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "none")


