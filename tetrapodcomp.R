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


#read data
Data <- read.csv("min20-2022.05.16.csv")
View(Data)



### Longevity model
#longevity neo
cutData <- Data[,c(5,9,10,11,13,40,42),drop=FALSE] 
cutData<-filter(cutData, Clade == "Mammalia")
cutData[cutData$max_longevity.months. < 0,] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevity.neo<-pglsSEyPagel(NeoplasiaPrevalence~max_longevity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevity.neo)

#grab r squared, lambda, and p values from summary 

r.v.longneo <- summary(longevity.neo)$corBeta
r.v.longneo <- format(r.v.longneo[2,1])
r.v.longneo <-signif(as.numeric(r.v.longneo)^2, digits= 2)
ld.v.longneo<- summary(longevity.neo)$modelStruct$corStruct
ld.v.longneo <- signif(ld.v.longneo[1], digits = 2)
p.v.longneo<-summary(longevity.neo)$tTable
p.v.longneo<-signif(p.v.longneo[2,4], digits = 3)

longneo<-ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(max_longevity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    limits = c(0,75),
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(longevity.neo)[1]*100, slope =  coef(longevity.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Max Longevity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  labs(title = "B")+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")


#ggsave(filename='longneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#ggsave(filename='wgtlong.pdf', width=9.5, height=18, limitsize=FALSE,bg="white")

#longevity mal
cutData <- Data[,c(5,9,10,11,17,40,42),drop=FALSE]
cutData<-filter(cutData, Clade == "Mammalia")
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevity.mal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity.months.),data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevity.mal)

#grab r squared, lambda, and p values from summary 

r.v.longmal <- summary(longevity.mal)$corBeta
r.v.longmal <- format(r.v.longmal[2,1])
r.v.longmal <-signif(as.numeric(r.v.longmal)^2, digits= 2)
ld.v.longmal<- summary(longevity.mal)$modelStruct$corStruct
ld.v.longmal <- signif(ld.v.longmal[1])
p.v.longmal<-summary(longevity.mal)$tTable
p.v.longmal<-signif(p.v.longmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(max_longevity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$max_longevity.months.)),log10(max(cutData$max_longevity.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(longevity.mal)[1]*100, slope =  coef(longevity.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("Max Longevity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Max Longevity",  
       subtitle =bquote(p-value:.(p.v.longmal)~R^2:.(r.v.longmal)~Lambda:.(ld.v.longmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=1.07, y=50.3, label = "5", size = 7)

#ggsave(filename='S5longmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")


### Longevity model. amph
#longevity neo

cutData <- Data[,c(5,9,10,11,13,40,42),drop=FALSE] 
cutData<-filter(cutData, Clade == "Amphibia")
cutData[cutData$max_longevity.months. < 0,] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevityamph.neo<-pglsSEyPagel(NeoplasiaPrevalence~max_longevity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevityamph.neo)

#grab r squared, lambda, and p values from summary 

r.v.longneoamph <- summary(longevityamph.neo)$corBeta
r.v.longneoamph <- format(r.v.longneoamph[2,1])
r.v.longneoamph <-signif(as.numeric(r.v.longneoamph)^2, digits= 2)
ld.v.longneoamph<- summary(longevityamph.neo)$modelStruct$corStruct
ld.v.longneoamph <- signif(ld.v.longneoamph[1], digits = 2)
p.v.longneoamph<-summary(longevityamph.neo)$tTable
p.v.longneoamph<-signif(p.v.longneoamph[2,4], digits = 3)



#longevity mal amph
cutData <- Data[,c(5,9,10,11,17,40,42),drop=FALSE] 
cutData<-filter(cutData, Clade == "Amphibia")
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevity.mal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity.months.),data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevity.mal)

#grab r squared, lambda, and p values from summary 

r.v.longmal <- summary(longevity.mal)$corBeta
r.v.longmal <- format(r.v.longmal[2,1])
r.v.longmal <-signif(as.numeric(r.v.longmal)^2, digits= 2)
ld.v.longmal<- summary(longevity.mal)$modelStruct$corStruct
ld.v.longmal <- signif(ld.v.longmal[1])
p.v.longmal<-summary(longevity.mal)$tTable
p.v.longmal<-signif(p.v.longmal[2,4], digits = 3)


### Longevity model. aves
#longevity neo
cutData <- Data[,c(4,5,9,10,11,13,40,42),drop=FALSE] 
cutData<-filter(cutData, Class == "Aves")
cutData[cutData$max_longevity.months. < 0,] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevityaves.neo<-pglsSEyPagel(NeoplasiaPrevalence~max_longevity.months.,data=cutData,
                                tree=pruned.tree,method="ML",se=SE)
summary(longevityaves.neo)

#grab r squared, lambda, and p values from summary 

r.v.longneoaves <- summary(longevityaves.neo)$corBeta
r.v.longneoaves <- format(r.v.longneoaves[2,1])
r.v.longneoaves <-signif(as.numeric(r.v.longneoaves)^2, digits= 2)
ld.v.longneoaves<- summary(longevityaves.neo)$modelStruct$corStruct
ld.v.longneoaves <- signif(ld.v.longneoaves[1], digits = 2)
p.v.longneoaves<-summary(longevityaves.neo)$tTable
p.v.longneoaves<-signif(p.v.longneoaves[2,4], digits = 3)



#longevity mal aves
cutData <- Data[,c(4,5,9,10,11,17,40,42),drop=FALSE] 
cutData<-filter(cutData, Class == "Aves")
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevityaves.mal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity.months.),data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevityaves.mal)

#grab r squared, lambda, and p values from summary 

r.v.longmalaves <- summary(longevityaves.mal)$corBeta
r.v.longmalaves <- format(r.v.longmalaves[2,1])
r.v.longmalaves <-signif(as.numeric(r.v.longmalaves)^2, digits= 2)
ld.v.longmalaves<- summary(longevityaves.mal)$modelStruct$corStruct
ld.v.longmalaves <- signif(ld.v.longmalaves[1])
p.v.longmalaves<-summary(longevityaves.mal)$tTable
p.v.longmalaves<-signif(p.v.longmalaves[2,4], digits = 3)




### Longevity model. squa
#longevity neo
cutData <- Data[,c(4,5,6,9,10,11,13,40,42),drop=FALSE] 
cutData<-filter(cutData, Orders == "Squamata")
cutData[cutData$max_longevity.months. < 0,] <-NA
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
longevitysqua.neo<-pglsSEyPagel(NeoplasiaPrevalence~max_longevity.months.,data=cutData,
                                tree=pruned.tree,method="ML",se=SE)
summary(longevitysqua.neo)

#grab r squared, lambda, and p values from summary 

r.v.longneosqua <- summary(longevitysqua.neo)$corBeta
r.v.longneosqua <- format(r.v.longneosqua[2,1])
r.v.longneosqua <-signif(as.numeric(r.v.longneosqua)^2, digits= 2)
ld.v.longneosqua<- summary(longevitysqua.neo)$modelStruct$corStruct
ld.v.longneosqua <- signif(ld.v.longneosqua[1], digits = 2)
p.v.longneosqua<-summary(longevitysqua.neo)$tTable
p.v.longneosqua<-signif(p.v.longneosqua[2,4], digits = 3)



#longevity mal squa
cutData <- Data[,c(4,5,6,9,10,11,17,40,42),drop=FALSE] 
cutData<-filter(cutData, Orders == "Squamata")
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

view(cutData)

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
longevitysqua.mal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity.months.),data=cutData,
                                tree=pruned.tree,method="ML",se=SE)
summary(longevitysqua.mal)

#grab r squared, lambda, and p values from summary 

r.v.longmalsqua <- summary(longevitysqua.mal)$corBeta
r.v.longmalsqua <- format(r.v.longmalsqua[2,1])
r.v.longmalsqua <-signif(as.numeric(r.v.longmalsqua)^2, digits= 2)
ld.v.longmalsqua<- summary(longevitysqua.mal)$modelStruct$corStruct
ld.v.longmalsqua <- signif(ld.v.longmalsqua[1])
p.v.longmalsqua<-summary(longevitysqua.mal)$tTable
p.v.longmalsqua<-signif(p.v.longmalsqua[2,4], digits = 3)


