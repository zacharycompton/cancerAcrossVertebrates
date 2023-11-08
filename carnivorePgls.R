library(nlme)
library(rms)
library(phytools)
library(geiger)
library(caper)
library(tidyverse)
library(cowplot)
library(ggrepel)
library(ggsci)
library(patchwork)

#make sure to run all of this before you get to work.
#pgls sey base (just run all of this)
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


#read data
Data <- read.csv("min20-2022.05.16.csv")
trophic<-read.csv("trophicData.csv")
View(Data)

Data<-left_join(Data, trophic, by = "Species")





#adult weight models
#adult weight neo

cutData <- Data[,c(5,6,9,10,11,13,38,42,43),drop=FALSE] 
cutData <- cutData %>%
  mutate(Carnivore = ifelse(WildTrophicLevel == "SecondaryCarnivore" |WildTrophicLevel == "PrimaryCarnivore", 1, 0))
cutData[cutData$adult_weight == -1, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
view(cutData)


#pgls model
adult.weight.neo.carn<-pglsSEyPagel(NeoplasiaPrevalence~log10(adult_weight.g.) + Carnivore,data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.neo.carn) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.neo.carn <- summary(adult.weight.neo.carn)$corBeta
r.v.adult.weight.neo.carn <- format(r.v.adult.weight.neo.carn[2,1])
r.v.adult.weight.neo.carn <-signif(as.numeric(r.v.adult.weight.neo.carn)^2, digits= 2)
ld.v.adult.weight.neo.carn<- summary(adult.weight.neo.carn)$modelStruct$corStruct
ld.v.adult.weight.neo.carn <- signif(ld.v.adult.weight.neo.carn[1], digits = 2)
p.v.adult.weight.neo.carn<-summary(adult.weight.neo.carn)$tTable
p.v.adult.weight.neo.carn<-signif(p.v.adult.weight.neo.carn[2,4], digits = 2)





#adult weight models
#adult weight mal

cutData <- Data[,c(5,6,9,10,11,17,38,42,43),drop=FALSE] 
cutData <- cutData %>%
  mutate(Carnivore = ifelse(WildTrophicLevel == "SecondaryCarnivore" |WildTrophicLevel == "PrimaryCarnivore", 1, 0))
cutData[cutData$adult_weight == -1, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
view(cutData)


#pgls model
adult.weight.mal.carn<-pglsSEyPagel(MalignancyPrevalence~log10(adult_weight.g.) + Carnivore,data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.mal.carn) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.mal.carn <- summary(adult.weight.mal.carn)$corBeta
r.v.adult.weight.mal.carn <- format(r.v.adult.weight.mal.carn[2,1])
r.v.adult.weight.mal.carn <-signif(as.numeric(r.v.adult.weight.mal.carn)^2, digits= 2)
ld.v.adult.weight.mal.carn<- summary(adult.weight.mal.carn)$modelStruct$corStruct
ld.v.adult.weight.mal.carn <- signif(ld.v.adult.weight.mal.carn[1], digits = 2)
p.v.adult.weight.mal.carn<-summary(adult.weight.mal.carn)$tTable
p.v.adult.weight.mal.carn<-signif(p.v.adult.weight.mal.carn[2,4], digits = 2)



#adult weight models
#adult weight neo

cutData <- Data[,c(5,6,9,10,11,13,30,38,42,43),drop=FALSE] 
cutData <- cutData %>%
  mutate(Carnivore = ifelse(WildTrophicLevel == "SecondaryCarnivore" |WildTrophicLevel == "PrimaryCarnivore", 1, 0))
cutData[cutData$adult_weight.g. == -1, ] <-NA
cutData <- na.omit(cutData)
cutData[cutData$Gestation.months. == -1, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
view(cutData)


#pgls model
adult.weight.gest.neo.carn<-pglsSEyPagel(NeoplasiaPrevalence~log10(adult_weight.g.) +log10(Gestation.months.)+ Carnivore,data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.gest.neo.carn) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.gest.neo.carn <- summary(adult.weight.gest.neo.carn)$corBeta
r.v.adult.weight.gest.neo.carn <- format(r.v.adult.weight.gest.neo.carn[2,1])
r.v.adult.weight.gest.neo.carn <-signif(as.numeric(r.v.adult.weight.gest.neo.carn)^2, digits= 2)
ld.v.adult.weight.gest.neo.carn<- summary(adult.weight.gest.neo.carn)$modelStruct$corStruct
ld.v.adult.weight.gest.neo.carn <- signif(ld.v.adult.weight.gest.neo.carn[1], digits = 2)
p.v.adult.weight.gest.neo.carn<-summary(adult.weight.gest.neo.carn)$tTable
p.v.adult.weight.gest.neo.carn<-signif(p.v.adult.weight.gest.neo.carn[2,4], digits = 2)

library(car)

covar<-vcov(adult.weight.gest.neo.carn)

# Assuming 'covar' is your covariance matrix
eigenvalues <- eigen(covar)$values

# Check for multicollinearity using eigenvalues
vif_values <- 1 / min(eigenvalues)

# Print or inspect VIF values
print(vif_values)







#adult weight models
#adult weight mal

cutData <- Data[,c(5,6,9,10,11,17,30,38,42,43),drop=FALSE] 
cutData <- cutData %>%
  mutate(Carnivore = ifelse(WildTrophicLevel == "SecondaryCarnivore" |WildTrophicLevel == "PrimaryCarnivore", 1, 0))
cutData[cutData$adult_weight == -1, ] <-NA
cutData <- na.omit(cutData)
cutData[cutData$Gestation.months. == -1, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
view(cutData)


#pgls model
adult.weight.gest.mal.carn<-pglsSEyPagel(MalignancyPrevalence~log10(adult_weight.g.) + log10(Gestation.months.)+Carnivore,data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.gest.mal.carn) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.mal.carn <- summary(adult.weight.mal.carn)$corBeta
r.v.adult.weight.mal.carn <- format(r.v.adult.weight.mal.carn[2,1])
r.v.adult.weight.mal.carn <-signif(as.numeric(r.v.adult.weight.mal.carn)^2, digits= 2)
ld.v.adult.weight.mal.carn<- summary(adult.weight.mal.carn)$modelStruct$corStruct
ld.v.adult.weight.mal.carn <- signif(ld.v.adult.weight.mal.carn[1], digits = 2)
p.v.adult.weight.mal.carn<-summary(adult.weight.mal.carn)$tTable
p.v.adult.weight.mal.carn<-signif(p.v.adult.weight.mal.carn[2,4], digits = 2)









