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
library(poolr)

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

foo <- function(x)
{
  nas <- is.na(x$coef)
  coef <- x$coef[!nas]
  cnames <- names(coef)
  coef <- matrix(rep(coef, 4), ncol = 4)
  dimnames(coef) <- list(cnames,
                         c("Estimate", "S.E.", "t", "Pr(T > |t|)"))
  df <- x$dfP - dim(coef)[1]
  coef[, 2] <- sqrt(diag(x$W))
  coef[, 3] <- coef[, 1]/coef[, 2]
  if (df < 0) {
    warning("not enough degrees of freedom to compute P-values.")
    coef[, 4] <- NA
  } else coef[, 4] <- 2 * (1 - pt(abs(coef[, 3]), df))
  coef
}


#read data
Data <- read.csv("min1LH.csv")
Data<-filter(Data, RecordsWithDenominators >= 1)
Data<-filter(Data, Class == "Mammalia")
View(Data)

#dataReport<-data.frame(min = NULL, model = NULL, pval = NULL, slope = NULL)

#adult weight models
#adult weight neo

cutData <- Data[,c(1,9,11,12,13,38),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min1.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
cutData$SE_simple <-1/sqrt(cutData$RecordsWithDenominators)
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]

#pgls model
adult.weight.neo<-pglsSEyPagel(NeoplasiaPrevalence~log(adult_weight),data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.neo) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.neo <- summary(adult.weight.neo)$corBeta
r.v.adult.weight.neo <- format(r.v.adult.weight.neo[2,1])
r.v.adult.weight.neo <-signif(as.numeric(r.v.adult.weight.neo)^2, digits= 2)
ld.v.adult.weight.neo<- summary(adult.weight.neo)$modelStruct$corStruct
ld.v.adult.weight.neo <- signif(ld.v.adult.weight.neo[1], digits = 2)
p.v.adult.weight.neo<-summary(adult.weight.neo)$tTable
p.v.adult.weight.neo<-signif(p.v.adult.weight.neo[2,4], digits = 2)
slope.adult.weight.neo<-summary(adult.weight.neo)$coefficients[2]



adultWeight<-log(cutData$adult_weight)
NeoplasiaOccurences<-cutData$NeoplasiaWithDenominators
NonOccurences<-(cutData$RecordsWithDenominators-cutData$NeoplasiaWithDenominators)

compar<-compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~adultWeight, phy = pruned.tree, family = "binomial")
compar


cNpval<-foo(compar)[2,4]
cNslope<-foo(compar)[2,1]

# ggplot(cutData, aes(x = adultWeight, y = NeoplasiaOccurences / (NeoplasiaOccurences + NonOccurences))) +
#   geom_point() +
#   labs(title = "Effect of Adult Weight on Neoplasia Occurrences", x = "Adult Weight", y = "Probability of Neoplasia Occurrence") +
#   theme_minimal()

#normal pgls
# Ensure 'species' column exists and is intended for matching
# Example assumes 'species' column in 'cutData' matches the tree tips in 'pruned.tree'
comp_data <- comparative.data(pruned.tree, cutData, "Species")

pgls<-pgls(NeoplasiaPrevalence~log(adult_weight), data = comp_data)



#adult weight mal
cutData <- Data[,c(5,9,10,11,17,38,27,16),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min1.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies <- cutData$Species

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
adult.weight.mal<-pglsSEyPagel(MalignancyPrevalence~log(adult_weight),data=cutData,
                               tree=pruned.tree,se=SE,method="ML")
summary(adult.weight.mal)

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.mal <- summary(adult.weight.mal)$corBeta
r.v.adult.weight.mal <- format(r.v.adult.weight.mal[2,1])
r.v.adult.weight.mal <-signif(as.numeric(r.v.adult.weight.mal)^2, digits= 2)
ld.v.adult.weight.mal<- summary(adult.weight.mal)$modelStruct$corStruct
ld.v.adult.weight.mal <- signif(ld.v.adult.weight.mal[1], digits= 2)
p.v.adult.weight.mal<-summary(adult.weight.mal)$tTable
p.v.adult.weight.mal<-signif(p.v.adult.weight.mal[2,4], digits = 3)
slope.adult.weight.mal<-summary(adult.weight.mal)$coefficients[2]



adultWeight<-log(cutData$adult_weight)
MalignantOccurences<-cutData$Malignant
NonOccurences<-(cutData$RecordsWithDenominators-cutData$Malignant)

comparMal<-compar.gee(cbind(MalignantOccurences, NonOccurences) ~adultWeight, phy = pruned.tree, family = "binomial")
comparMal
cMpval<-foo(comparMal)[2,4]
cMslope<-foo(comparMal)[2,1]


# ggplot(cutData, aes(x = adultWeight, y = NeoplasiaOccurences / (NeoplasiaOccurences + NonOccurences))) +
#   geom_point() +
#   labs(title = "Effect of Adult Weight on Neoplasia Occurrences", x = "Adult Weight", y = "Probability of Neoplasia Occurrence") +
#   theme_minimal()

#normal pgls
# Ensure 'species' column exists and is intended for matching
# Example assumes 'species' column in 'cutData' matches the tree tips in 'pruned.tree'
comp_dataMal <- comparative.data(pruned.tree, cutData, "Species")

pglsMal<-pgls(MalignancyPrevalence~log(adult_weight), data = comp_dataMal)



pglsPagel_row<-data.frame(min = c(min(Data$RecordsWithDenominators),min(Data$RecordsWithDenominators)), model = c("PGLSSEYneo","PGLSSEYmal"), pval = c(p.v.adult.weight.neo, p.v.adult.weight.mal), slope = c(slope.adult.weight.neo,slope.adult.weight.mal))



compar_row<-data.frame(min = c(min(Data$RecordsWithDenominators),min(Data$RecordsWithDenominators)), model = c("Comparneo","Comparmal"), pval = c(cNpval, cMpval), slope = c(cNslope,cMslope))


dataReport <- rbind(dataReport, pglsPagel_row)
dataReport <- rbind(dataReport, compar_row)

#long
cutData <- Data[,c(1,9,11,12,13,38,40,27),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min1.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]

max_longevity.months.<-log(cutData$max_longevity)
NeoplasiaOccurences<-cutData$NeoplasiaWithDenominators
NonOccurences<-(cutData$RecordsWithDenominators-cutData$NeoplasiaWithDenominators)

compar<-compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~max_longevity.months., phy = pruned.tree, family = "binomial")
compar


#long weight mal
cutData <- Data[,c(5,9,16,17,38,27,16,40,11),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min1.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies <- cutData$Species

pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
view(cutData)


max_longevity.months.<-log(cutData$max_longevity)
MalignantOccurences<-cutData$Malignant
NonOccurences<-(cutData$RecordsWithDenominators-cutData$Malignant)

comparMal<-compar.gee(cbind(MalignantOccurences, NonOccurences) ~max_longevity.months., phy = pruned.tree, family = "binomial")
comparMal
cMpval<-foo(comparMal)[2,4]
cMslope<-foo(comparMal)[2,1]




#gest
cutData <- Data[,c(1,9,11,12,13,38,30,27),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min1.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]

Gestation<-log(cutData$Gestation)
NeoplasiaOccurences<-cutData$NeoplasiaWithDenominators
NonOccurences<-(cutData$RecordsWithDenominators-cutData$NeoplasiaWithDenominators)

compar<-compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~Gestation, phy = pruned.tree, family = "binomial")
compar


#gest weight mal
cutData <- Data[,c(5,9,16,17,38,27,16,30,11),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min1.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies <- cutData$Species

pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simple,cutData$Species)[rownames(cutData)]
view(cutData)


Gestation<-log(cutData$Gestation)
MalignantOccurences<-cutData$Malignant
NonOccurences<-(cutData$RecordsWithDenominators-cutData$Malignant)

comparMal<-compar.gee(cbind(MalignantOccurences, NonOccurences) ~Gestation, phy = pruned.tree, family = "binomial")
comparMal
cMpval<-foo(comparMal)[2,4]
cMslope<-foo(comparMal)[2,1]
