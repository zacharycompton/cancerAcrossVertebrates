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
View(Data)




#adult weight models
#adult weight neo

cutData <- Data[,c(5,9,10,11,13,38,42),drop=FALSE] 
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
adult.weight.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(adult_weight.g.),data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.neo) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.neo <- summary(adult.weight.neo)$corBeta
r.v.adult.weight.neo <- format(r.v.adult.weight.neo[2,1])
r.v.adult.weight.neo <-signif(as.numeric(r.v.adult.weight.neo)^2, digits= 2)
ld.v.adult.weight.neo<- summary(adult.weight.neo)$modelStruct$corStruct
ld.v.adult.weight.neo <- signif(ld.v.adult.weight.neo[1], digits = 2)
p.v.adult.weight.neo<-summary(adult.weight.neo)$tTable
p.v.adult.weight.neo<-signif(p.v.adult.weight.neo[2,4], digits = 2)


#plot
wgtneo<-ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(adult_weight.g.)))+
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ),)+
  scale_y_continuous(
    limits = c(0,75),
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(adult.weight.neo)[1]*100, slope =  coef(adult.weight.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Adult Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  geom_text_repel(aes(label=ifelse( NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(title = "A")+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")


#ggsave(filename='wgtneol.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#adult weight mal
cutData <- Data[,c(5,9,10,11,17,38,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
adult.weight.mal<-pglsSEyPagel(MalignancyPrevalence~log10(adult_weight.g.),data=cutData,
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

#remove limits from scale y continous, add last two lines with label
ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(adult_weight.g.)))+
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$adult_weight.g.)),log10(max(cutData$adult_weight.g.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(adult.weight.mal)[1]*100, slope =  coef(adult.weight.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Adult Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  labs(title = "Malignancy Prevalence vs. Adult Weight",
       subtitle =bquote(p-value:.(p.v.adult.weight.mal)~R^2:.(r.v.adult.weight.mal)~Lambda:.(ld.v.adult.weight.mal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=.42, y=50.3, label = "1", size = 7)



#ggsave(filename='S1wgtmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#gestation models
#gestation neo
cutData <- Data[,c(5,9,10,11,13,30,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
gestation.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(Gestation.months.),data=cutData,
                            tree=pruned.tree,method="ML",se=SE)

summary(gestation.neo)

#grab r squared, lambda, and p values from summary 

r.v.gestneo <- summary(gestation.neo)$corBeta
r.v.gestneo <- format(r.v.gestneo[2,1])
r.v.gestneo<-signif(as.numeric(r.v.gestneo)^2, digits= 2)
ld.v.gestneo<- summary(gestation.neo)$modelStruct$corStruct
ld.v.gestneo <- signif(ld.v.gestneo[1], digits = 2)
p.v.gestneo<-summary(gestation.neo)$tTable
p.v.gestneo<-signif(p.v.gestneo[2,4], digits = 3)

gestneo<-ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(Gestation.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    limits = c(0,75),
    breaks = c(0, 25,45,75),
    labels = c(0, 25,45,75))+
  geom_abline(intercept = coef(gestation.neo)[1]*100, slope =  coef(gestation.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Gestation") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  labs(title="A")

##ggsave(filename='S2gestneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#gestation mal
cutData <- Data[,c(5,9,10,11,17,30,42),drop=FALSE] 
cutData[cutData$Gestation.months. < 0, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
gestation.mal<-pglsSEyPagel(MalignancyPrevalence~log10(Gestation.months.),data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(gestation.mal)

#grab r squared, lambda, and p values from summary 

r.v.gestmal <- summary(gestation.mal)$corBeta
r.v.gestmal <- format(r.v.gestmal[2,1])
r.v.gestmal<-signif(as.numeric(r.v.gestmal)^2, digits= 2)
ld.v.gestmal<- summary(gestation.mal)$modelStruct$corStruct
ld.v.gestmal<- signif(ld.v.gestmal[1], digits = 2)
p.v.gestmal<-summary(gestation.mal)$tTable
p.v.gestmal<-signif(p.v.gestmal[2,4], digits = 3)

gestmal<-ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(Gestation.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    limits = c(0,75),
    breaks = c(0, 25,45,75),
    labels = c(0, 25,45,75))+
  geom_abline(intercept = coef(gestation.mal)[1]*100, slope =  coef(gestation.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Gestation") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  labs(title="B")

##ggsave(filename='gestmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

gestneo/gestmal

#ggsave(filename='gestneomal.pdf', width=9.5, height=18, limitsize=FALSE,bg="white")


#litter size models 
#litter size neo
cutData <- Data[,c(5,9,10,11,13,33,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
litter.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(litter_size),data=cutData,
                         tree=pruned.tree,method="ML",se=SE)
summary(litter.neo)

#grab r squared, lambda, and p values from summary 

r.v.litneo <- summary(litter.neo)$corBeta
r.v.litneo <- format(r.v.litneo[2,1])
r.v.litneo <-signif(as.numeric(r.v.litneo)^2, digits= 2)
ld.v.litneo<- summary(litter.neo)$modelStruct$corStruct
ld.v.litneo <- signif(ld.v.litneo[1], digits = 2)
p.v.litneo<-summary(litter.neo)$tTable
p.v.litneo<-signif(p.v.litneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(litter_size))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(litter.neo)[1]*100, slope =  coef(litter.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Litter Size") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Litter Size",  
       subtitle =bquote(p-value:.(p.v.litneo)~R^2:.(r.v.litneo)~Lambda:.(ld.v.litneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$litter_size)),log10(max(cutData$litter_size))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=-0.35, y=83.8, label = "3", size = 7)


#ggsave(filename='S3litneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#litter size mal
cutData <- Data[,c(5,9,10,11,17,33,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)

tree <- read.tree("min20Fixed516.nwk")
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
litter.mal <- pglsSEyPagel(MalignancyPrevalence~log10(litter_size),data=cutData,
                           tree=pruned.tree,method="ML",se=SE)
summary(litter.mal)

#grab r squared, lambda, and p values from summary 

r.v.litmal <- summary(litter.mal)$corBeta
r.v.litmal <- format(r.v.litmal[2,1])
r.v.litmal <-signif(as.numeric(r.v.litmal)^2, digits= 2)
ld.v.litmal<- summary(litter.mal)$modelStruct$corStruct
ld.v.litmal <- signif(ld.v.litmal[1], digits = 2)
p.v.litmal<-summary(litter.mal)$tTable
p.v.litmal<-signif(p.v.litmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(litter_size))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$litter_size)),log10(max(cutData$litter_size))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(litter.mal)[1]*100, slope =  coef(litter.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Litter Size") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Litter Size",  
       subtitle =bquote(p-value:.(p.v.litmal)~R^2:.(r.v.litmal)~Lambda:.(ld.v.litmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=-0.35, y=50.3, label = "4", size = 7)

#ggsave(filename='S4litmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")


### Longevity model
#longevity neo
cutData <- Data[,c(5,9,10,11,13,40,42),drop=FALSE] 
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


#create weight over longevity model
wgtneo/longneo
#ggsave(filename='wgtlong.pdf', width=9.5, height=18, limitsize=FALSE,bg="white")

#longevity mal

cutData <- Data[,c(5,9,10,11,17,40,42),drop=FALSE] 
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


##BMR models
#bmr neo
cutData <- Data[,c(5,9,10,11,13,41,42),drop=FALSE] 
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
BMR.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(metabolic_rate),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)
summary(BMR.neo)

#grab r squared, lambda, and p values from summary 

r.v.bmrneo <- summary(BMR.neo)$corBeta
r.v.bmrneo <- format(r.v.bmrneo[2,1])
r.v.bmrneo <-signif(as.numeric(r.v.bmrneo)^2, digits= 2)
ld.v.bmrneo<- summary(BMR.neo)$modelStruct$corStruct
ld.v.bmrneo <- signif(ld.v.bmrneo[1], digits = 2)
p.v.bmrneo<-summary(BMR.neo)$tTable
p.v.bmrneo<-signif(p.v.bmrneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(metabolic_rate))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(BMR.neo)[1]*100, slope =  coef(BMR.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Basal Metabolic Rate") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Metabolic Rate in Mammals",  
       subtitle =bquote(p-value:.(p.v.bmrneo)~R^2:.(r.v.bmrneo)~Lambda:.(ld.v.bmrneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$metabolic_rate)),log10(max(cutData$metabolic_rate))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=.9, y=83.8, label = "6", size = 7)

#ggsave(filename='S6bmrneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#bmr mal

cutData <- Data[,c(5,9,10,11,17,41,42),drop=FALSE] 
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
BMR.mal<-pglsSEyPagel(MalignancyPrevalence~log10(metabolic_rate),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)
summary(BMR.mal)

#grab r squared, lambda, and p values from summary 

r.v.bmrmal <- summary(BMR.mal)$corBeta
r.v.bmrmal <- format(r.v.bmrmal[2,1])
r.v.bmrmal <-signif(as.numeric(r.v.bmrmal)^2, digits= 2)
ld.v.bmrmal<- summary(BMR.mal)$modelStruct$corStruct
ld.v.bmrmal <- signif(ld.v.bmrmal[1], digits = 2)
p.v.bmrmal<-summary(BMR.mal)$tTable
p.v.bmrmal<-signif(p.v.bmrmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(metabolic_rate))) + 
  geom_abline(intercept = coef(BMR.mal)[1]*100, slope =  coef(BMR.mal)[2]*100,
              color = 'grey',size = 1.2) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$metabolic_rate)),log10(max(cutData$metabolic_rate))),
                  ylim = c(0,45),clip = "off")+
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Basal Metabolic Rate") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Metabolic Rate",  
       subtitle =bquote(p-value:.(p.v.bmrmal)~R^2:.(r.v.bmrmal)~Lambda:.(ld.v.bmrmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=.9, y=50.3, label = "7", size = 7)

#ggsave(filename='S7bmrmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#wxl models #weight and longevity
#wxl neo
cutData <- Data[,c(5,9,10,11,13,40,38,42),drop=FALSE] 
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
wxl.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(max_longevity.months.*adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wxl.neo)

#grab r squared, lambda, and p values from summary 

r.v.wxneo <- summary(wxl.neo)$corBeta
r.v.wxneo <- format(r.v.wxneo[2,1])
r.v.wxneo <-signif(as.numeric(r.v.wxneo)^2, digits= 2)
ld.v.wxneo<- summary(wxl.neo)$modelStruct$corStruct
ld.v.wxneo <- signif(ld.v.wxneo[1], digits = 2)
p.v.wxneo<-summary(wxl.neo)$tTable
p.v.wxneo<-signif(p.v.wxneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(adult_weight.g.*max_longevity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(wxl.neo)[1]*100, slope =  coef(wxl.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Adult Weight(g)*Longevity(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Max Longevity*Weight",  
       subtitle =bquote(p-value:.(p.v.wxneo)~R^2:.(r.v.wxneo)~Lambda:.(ld.v.wxneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.*cutData$max_longevity.months.))),max(log10(cutData$adult_weight.g.*cutData$max_longevity.months.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=2, y=83.8, label = "8", size = 7)

#ggsave(filename='S8wgtlongneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#wxl mal
cutData <- Data[,c(5,9,10,11,17,40,38,42),drop=FALSE] 
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
wxl.mal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity.months.*adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wxl.mal)

#grab r squared, lambda, and p values from summary 

r.v.wxmal <- summary(wxl.mal)$corBeta
r.v.wxmal <- format(r.v.wxmal[2,1])
r.v.wxmal <-signif(as.numeric(r.v.wxmal)^2, digits= 2)
ld.v.wxmal<- summary(wxl.mal)$modelStruct$corStruct
ld.v.wxmal <- signif(ld.v.wxmal[1], digits = 2)
p.v.wxmal<-summary(wxl.mal)$tTable
p.v.wxmal<-signif(p.v.wxmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(adult_weight.g.*max_longevity.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.*cutData$max_longevity.months.))),max(log10(cutData$adult_weight.g.*cutData$max_longevity.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(wxl.mal)[1]*100, slope =  coef(wxl.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Adult Weight(g)*Longevity(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Max Longevity*Weight",  
       subtitle =bquote(p-value:.(p.v.wxmal)~R^2:.(r.v.wxmal)~Lambda:.(ld.v.wxmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  annotate("text", x=2, y=50.3, label = "9", size = 7)

#ggsave(filename='S9wgtlongmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#litters per year neo
cutData <- Data[,c(5,9,10,11,13,34,42),drop=FALSE] 
cutData <- cutData[!(cutData$common_name %in% c("Dwarf caimen","Agassiz's desert tortoise","Egyptian snouted cobra","Mangrove snake","Tiger rat snake","Madagscar tree boa")),]
cutData[cutData$litters_year < 0, ] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
lityear.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(litters_year),data=cutData,
                          tree=pruned.tree,method="ML",se=SE)
summary(lityear.neo)

#grab r squared, lambda, and p values from summary 

r.v.lyearneo <- summary(lityear.neo)$corBeta
r.v.lyearneo <- format(r.v.lyearneo[2,1])
r.v.lyearneo<-signif(as.numeric(r.v.lyearneo)^2, digits= 2)
ld.v.lyearneo<- summary(lityear.neo)$modelStruct$corStruct
ld.v.lyearneo <- signif(ld.v.lyearneo[1], digits = 2)
p.v.lyearneo<-summary(lityear.neo)$tTable
p.v.lyearneo<-signif(p.v.lyearneo[2,4], digits = 3)

lityearneo<-ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(litters_year))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    limits = c(0,75),
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(lityear.neo)[1]*100, slope =  coef(lityear.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Litters per Year") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5)+
  labs(title = "A")+
  guides(colour = guide_legend(override.aes=list(size=5),order = 1) ,size= guide_legend(order=2)) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(size="Total Necropsies",colour="Clade")

#ggsave(filename='lityearneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#create lit year over gest model
lityearneo/gestmal
#ggsave(filename='litgest.pdf', width=9.5, height=18, limitsize=FALSE,bg="white")

#litter per year mal
cutData <- Data[,c(5,9,10,11,17,34,42),drop=FALSE] 
cutData <- cutData[!(cutData$common_name %in% c("Dwarf caimen","Agassiz's desert tortoise","Egyptian snouted cobra","Mangrove snake","Tiger rat snake","Madagscar tree boa")),]
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
lityear.mal<-pglsSEyPagel(MalignancyPrevalence~log10(litters_year),data=cutData,
                          tree=pruned.tree,method="ML",se=SE)
summary(lityear.mal)

#grab r squared, lambda, and p values from summary 

r.v.lyearmal <- summary(lityear.mal)$corBeta
r.v.lyearmal <- format(r.v.lyearmal[2,1])
r.v.lyearmal<-signif(as.numeric(r.v.lyearmal)^2, digits= 2)
ld.v.lyearmal<- summary(lityear.mal)$modelStruct$corStruct
ld.v.lyearmal<- signif(ld.v.lyearmal[1], digits = 2)
p.v.lyearmal<-summary(lityear.mal)$tTable
p.v.lyearmal<-signif(p.v.lyearmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(litters_year))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$litters_year)),log10(max(cutData$litters_year))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(lityear.mal)[1]*100, slope =  coef(lityear.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Litters per Year") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Litters per Year",  
       subtitle =bquote(p-value:.(p.v.lyearmal)~R^2:.(r.v.lyearmal)~Lambda:.(ld.v.lyearmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text",  x=-.85, y=50.3, label = "10", size = 7)

#ggsave(filename='S10lityearmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

####Female Maturity neo
cutData <- Data[,c(5,9,10,11,13,28,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
Fmaturity.neo<-pglsSEyPagel(NeoplasiaPrevalence~female_maturity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(Fmaturity.neo)

#grab r squared, lambda, and p values from summary 

r.v.fmaturityneo <- summary(Fmaturity.neo)$corBeta
r.v.fmaturityneo <- format(r.v.fmaturityneo[2,1])
r.v.fmaturityneo<-signif(as.numeric(r.v.fmaturityneo)^2, digits= 2)
ld.v.fmaturityneo<- summary(Fmaturity.neo)$modelStruct$corStruct
ld.v.fmaturityneo <- signif(ld.v.fmaturityneo[1], digits = 2)
p.v.fmaturityneo<-summary(Fmaturity.neo)$tTable
p.v.fmaturityneo<-signif(p.v.fmaturityneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(female_maturity.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(Fmaturity.neo)[1]*100, slope =  coef(Fmaturity.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Female Maturity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Female Maturity",  
       subtitle =bquote(p-value:.(p.v.fmaturityneo)~R^2:.(r.v.fmaturityneo)~Lambda:.(ld.v.fmaturityneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$female_maturity.months.)),log10(max(cutData$female_maturity.months.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=-.09, y=83.9, label = "11", size = 7)

#ggsave(filename='S11femmatneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#Female Maturity Mal
cutData <- Data[,c(5,9,10,11,17,28,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
Fmaturity.mal<-pglsSEyPagel(MalignancyPrevalence~female_maturity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(Fmaturity.mal)

#grab r squared, lambda, and p values from summary 

r.v.fmaturitymal <- summary(Fmaturity.mal)$corBeta
r.v.fmaturitymal <- format(r.v.fmaturitymal[2,1])
r.v.fmaturitymal<-signif(as.numeric(r.v.fmaturitymal)^2, digits= 2)
ld.v.fmaturitymal<- summary(Fmaturity.mal)$modelStruct$corStruct
ld.v.fmaturitymal <- signif(ld.v.fmaturitymal[1], digits = 2)
p.v.fmaturitymal<-summary(Fmaturity.mal)$tTable
p.v.fmaturitymal<-signif(p.v.fmaturitymal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(female_maturity.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$female_maturity.months.)),log10(max(cutData$female_maturity.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(Fmaturity.mal)[1]*100, slope =  coef(Fmaturity.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("Female Maturity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Female Maturity",  
       subtitle =bquote(p-value:.(p.v.fmaturitymal)~R^2:.(r.v.fmaturitymal)~Lambda:.(ld.v.fmaturitymal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=-.09, y=50.3, label = "12", size = 7)

#ggsave(filename='S12femmatmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

####Male Maturity neo
cutData <- Data[,c(5,9,10,11,13,29,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
Mmaturity.neo<-pglsSEyPagel(NeoplasiaPrevalence~male_maturity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(Mmaturity.neo)

#grab r squared, lambda, and p values from summary 

r.v.Mmaturityneo <- summary(Mmaturity.neo)$corBeta
r.v.Mmaturityneo <- format(r.v.Mmaturityneo[2,1])
r.v.Mmaturityneo<-signif(as.numeric(r.v.Mmaturityneo)^2, digits= 2)
ld.v.Mmaturityneo<- summary(Mmaturity.neo)$modelStruct$corStruct
ld.v.Mmaturityneo <- signif(ld.v.Mmaturityneo[1], digits = 2)
p.v.Mmaturityneo<-summary(Mmaturity.neo)$tTable
p.v.Mmaturityneo<-signif(p.v.Mmaturityneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(male_maturity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(Mmaturity.neo)[1]*100, slope =  coef(Mmaturity.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("Male Maturity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Male Maturity",  
       subtitle =bquote(p-value:.(p.v.Mmaturityneo)~R^2:.(r.v.Mmaturityneo)~Lambda:.(ld.v.Mmaturityneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$male_maturity.months.)),log10(max(cutData$male_maturity.months.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=-0.74, y=83.8, label = "13", size = 7)

#ggsave(filename='S13malematneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#Male Maturity Mal
cutData <- Data[,c(5,9,10,11,17,29,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
Mmaturity.mal<-pglsSEyPagel(MalignancyPrevalence~male_maturity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(Mmaturity.mal)

#grab r squared, lambda, and p values from summary 

r.v.Mmaturitymal <- summary(Mmaturity.mal)$corBeta
r.v.Mmaturitymal <- format(r.v.Mmaturitymal[2,1])
r.v.Mmaturitymal<-signif(as.numeric(r.v.Mmaturitymal)^2, digits= 2)
ld.v.Mmaturitymal<- summary(Mmaturity.mal)$modelStruct$corStruct
ld.v.Mmaturitymal <- signif(ld.v.Mmaturitymal[1], digits = 2)
p.v.Mmaturitymal<-summary(Mmaturity.mal)$tTable
p.v.Mmaturitymal<-signif(p.v.Mmaturitymal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(male_maturity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$male_maturity.months.)),log10(max(cutData$male_maturity.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(Fmaturity.mal)[1]*100, slope =  coef(Fmaturity.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("Male Maturity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Male Maturity",  
       subtitle =bquote(p-value:.(p.v.Mmaturitymal)~R^2:.(r.v.Mmaturitymal)~Lambda:.(ld.v.Mmaturitymal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=-.74, y=50.3, label = "14", size = 7)

#ggsave(filename='S14malematmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#Weaning Weight neo
cutData <- Data[,c(5,9,10,11,13,37,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
weanw.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(weaning_weight.g.),data=cutData,
                        tree=pruned.tree,method="ML",se=SE)
summary(weanw.neo)

#grab r squared, lambda, and p values from summary 

r.v.weanwneo <- summary(weanw.neo)$corBeta
r.v.weanwneo <- format(r.v.weanwneo[2,1])
r.v.weanwneo<-signif(as.numeric(r.v.weanwneo)^2, digits= 2)
ld.v.weanwneo<- summary(weanw.neo)$modelStruct$corStruct
ld.v.weanwneo <- signif(ld.v.weanwneo[1], digits = 2)
p.v.weanwneo<-summary(weanw.neo)$tTable
p.v.weanwneo<-signif(p.v.weanwneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(weaning_weight.g.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(weanw.neo)[1]*100, slope =  coef(weanw.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Weaning Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Weaning Weight",  
       subtitle =bquote(p-value:.(p.v.weanwneo)~R^2:.(r.v.weanwneo)~Lambda:.(ld.v.weanwneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$weaning_weight.g.)),log10(max(cutData$weaning_weight.g.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text",  x=-.64, y=83.8, label = "15", size = 7)

#ggsave(filename='S15weanneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#Weaning Weight Mal
cutData <- Data[,c(5,9,10,11,17,37,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
weanw.mal<-pglsSEyPagel(MalignancyPrevalence~log10(weaning_weight.g.),data=cutData,
                        tree=pruned.tree,method="ML",se=SE)
summary(weanw.mal)

#grab r squared, lambda, and p values from summary 

r.v.weanwmal <- summary(weanw.mal)$corBeta
r.v.weanwmal <- format(r.v.weanwmal[2,1])
r.v.weanwmal<-signif(as.numeric(r.v.weanwmal)^2, digits= 2)
ld.v.weanwmal<- summary(weanw.mal)$modelStruct$corStruct
ld.v.weanwmal <- signif(ld.v.weanwmal[1], digits = 2)
p.v.weanwmal<-summary(weanw.mal)$tTable
p.v.weanwmal<-signif(p.v.weanwmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(weaning_weight.g.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$weaning_weight.g.)),log10(max(cutData$weaning_weight.g.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(Fmaturity.mal)[1]*100, slope =  coef(Fmaturity.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Weaning Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Weaning Weight",  
       subtitle =bquote(p-value:.(p.v.weanwmal)~R^2:.(r.v.weanwmal)~Lambda:.(ld.v.weanwmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=-.64, y=50.3, label = "16", size = 7)

#ggsave(filename='S16weanmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#Growth Rate Neo
cutData <- Data[,c(5,9,10,11,13,39,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
GrowthR.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(growth_rate.1.days.),data=cutData,
                          tree=pruned.tree,method="ML",se=SE)
summary(GrowthR.neo)

#grab r squared, lambda, and p values from summary 

r.v.GrowthRneo <- summary(GrowthR.neo)$corBeta
r.v.GrowthRneo <- format(r.v.GrowthRneo[2,1])
r.v.GrowthRneo<-signif(as.numeric(r.v.GrowthRneo)^2, digits= 2)
ld.v.GrowthRneo<- summary(GrowthR.neo)$modelStruct$corStruct
ld.v.GrowthRneo <- signif(ld.v.GrowthRneo[1], digits = 2)
p.v.GrowthRneo<-summary(GrowthR.neo)$tTable
p.v.GrowthRneo<-signif(p.v.GrowthRneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(growth_rate.1.days.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(GrowthR.neo)[1]*100, slope =  coef(GrowthR.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Growth Rate") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Growth Rate",  
       subtitle =bquote(p-value:.(p.v.GrowthRneo)~R^2:.(r.v.GrowthRneo)~Lambda:.(ld.v.GrowthRneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$growth_rate.1.days.)),log10(max(cutData$growth_rate.1.days.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text",x=-3.36, y=83.8, label = "17", size = 7)

#ggsave(filename='S17growneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#Growth Rate Mal
cutData <- Data[,c(5,9,10,11,17,39,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

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
GrowthR.mal<-pglsSEyPagel(MalignancyPrevalence~log10(growth_rate.1.days.),data=cutData,
                          tree=pruned.tree,method="ML",se=SE)
summary(GrowthR.mal)

#grab r squared, lambda, and p values from summary 

r.v.GrowthRmal <- summary(GrowthR.mal)$corBeta
r.v.GrowthRmal <- format(r.v.GrowthRmal[2,1])
r.v.GrowthRmal<-signif(as.numeric(r.v.GrowthRmal)^2, digits= 2)
ld.v.GrowthRmal<- summary(GrowthR.mal)$modelStruct$corStruct
ld.v.GrowthRmal <- signif(ld.v.GrowthRmal[1], digits = 2)
p.v.GrowthRmal<-summary(GrowthR.mal)$tTable
p.v.GrowthRmal<-signif(p.v.GrowthRmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(growth_rate.1.days.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c(log10(min(cutData$growth_rate.1.days.)),log10(max(cutData$growth_rate.1.days.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(GrowthR.mal)[1]*100, slope =  coef(GrowthR.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Growth Rate") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Growth Rate",  
       subtitle =bquote(p-value:.(p.v.GrowthRmal)~R^2:.(r.v.GrowthRmal)~Lambda:.(ld.v.GrowthRmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  annotate("text", x=-3.35, y=50.3, label = "18", size = 7)

#ggsave(filename='S18growmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#w+g

cutData <- Data[,c(5,9,10,11,13,38,30,42),drop=FALSE] 
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
wpl.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(Gestation.months.)+log10(adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

coef(wpl.neo)

summary(wpl.neo)

#grab r squared, lambda, and p values from summary 

r.v.wpneo <- summary(wpl.neo)$corBeta
r.v.wpneo <- format(r.v.wpneo[3,1])
r.v.wpneo <-signif(as.numeric(r.v.wpneo)^2, digits= 2)
ld.v.wpneo<- summary(wpl.neo)$modelStruct$corStruct
ld.v.wpneo <- signif(ld.v.wpneo[1], digits = 2)
p.v.wpneo<-summary(wpl.neo)$tTable
p.v.wpneogest<-signif(p.v.wpneo[2,4], digits = 3)
p.v.wpneowgt<-signif(p.v.wpneo[3,4], digits = 2)
p.v.wpgneogest<-signif(p.v.wpneo[2,4], digits = 3)
p.v.wpgneolong<-signif(p.v.wpneo[3,4], digits = 2)
c.wpplneogest<-coef(wpl.neo)[2]
c.wpplneolong<-coef(wpl.neo)[3]



pvalues<-c(p.v.wpgneogest,p.v.wpgneolong)

combopwpgneo<-fisher(pvalues)


ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(adult_weight.g.)+log10(Gestation.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(wpl.neo)[1]*100, slope =  coef(wpl.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Adult Weight(g)+Gestation(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Max Longevity+Gestation",  
       subtitle =bquote(p-value:.(p.v.wpneo)~R^2:.(r.v.wpneo)~Lambda:.(ld.v.wpneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=.13, y=83.8, label = "19", size = 7)

#ggsave(filename='S19wgtgestneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")



#w+g mal

cutData <- Data[,c(5,9,10,11,17,38,30,42),drop=FALSE] 
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
wpl.mal<-pglsSEyPagel(MalignancyPrevalence~log10(Gestation.months.)+log10(adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.mal)

#grab r squared, lambda, and p values from summary 

r.v.wpmal <- summary(wpl.mal)$corBeta
r.v.wpmal <- format(r.v.wpmal[2,1])
r.v.wpmal <-signif(as.numeric(r.v.wpmal)^2, digits= 2)
ld.v.wpmal<- summary(wpl.mal)$modelStruct$corStruct
ld.v.wpmal <- signif(ld.v.wpmal[1], digits = 2)
p.v.wpmal<-summary(wpl.mal)$tTable
p.v.wpgmalgest<-signif(p.v.wpmal[2,4], digits = 3)
p.v.wpgmallong<-signif(p.v.wpmal[3,4], digits = 2)
c.wpplmalgest<-coef(wpl.mal)[2]
c.wpplmallong<-coef(wpl.mal)[3]



pvalues<-c(p.v.wpgmalgest,p.v.wpgmallong)

combopwpgmal<-fisher(pvalues)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(Gestation.months.)+log10(adult_weight.g.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(wpl.mal)[1]*100, slope =  coef(wpl.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Adult Weight(g)+Gestation(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Max Longevity+Gestation",  
       subtitle =bquote(p-value:.(p.v.wpmal)~R^2:.(r.v.wpmal)~Lambda:.(ld.v.wpmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")+
  annotate("text", x=.13, y=50.3, label = "20", size = 7)


#ggsave(filename='S20wgtgestmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#g+l

cutData <- Data[,c(5,9,10,11,13,40,30,42),drop=FALSE] 
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
wppl.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(Gestation.months.)+log10(max_longevity.months.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wppl.neo)
coef(wppl.neo)

#grab r squared, lambda, and p values from summary 

r.v.wpplneo <- summary(wppl.neo)$corBeta
r.v.wpplneo <- format(r.v.wpplneo[2,1])
r.v.wpplneo <-signif(as.numeric(r.v.wpplneo)^2, digits= 2)
ld.v.wpplneo<- summary(wppl.neo)$modelStruct$corStruct
ld.v.wpplneo <- signif(ld.v.wpplneo[1], digits = 2)
p.v.wpplneo<-summary(wppl.neo)$tTable
p.v.wpplneogest<-signif(p.v.wpplneo[2,4], digits = 3)
p.v.wpplneolong<-signif(p.v.wpplneo[3,4], digits = 2)
c.wpplneogest<-coef(wppl.neo)[2]
c.wpplneolong<-coef(wppl.neo)[3]


#combine p values for fisher p value
pvalues<-c(p.v.wpplneogest,p.v.wpplneolong)

combopwpplneo<-fisher(pvalues)



#g+l mal

cutData <- Data[,c(5,9,10,11,17,40,30,42),drop=FALSE] 
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
wppl.mal<-pglsSEyPagel(MalignancyPrevalence~log10(Gestation.months.)+log10(max_longevity.months.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wppl.mal)

#grab r squared, lambda, and p values from summary 

r.v.wpplmal <- summary(wppl.mal)$corBeta
r.v.wpplmal <- format(r.v.wpplmal[2,1])
r.v.wpplmal <-signif(as.numeric(r.v.wpplmal)^2, digits= 2)
ld.v.wpplmal<- summary(wppl.mal)$modelStruct$corStruct
ld.v.wpplmal <- signif(ld.v.wpplmal[1], digits = 2)
p.v.wpplmal<-summary(wppl.mal)$tTable
p.v.wpplmalgest<-signif(p.v.wpplmal[2,4], digits = 3)
p.v.wpplmallong<-signif(p.v.wpplmal[3,4], digits = 2)
c.wpplmalgest<-coef(wppl.mal)[2]
c.wpplmallong<-coef(wppl.mal)[3]


#combine p values for fisher p value
pvalues<-c(p.v.wpplmalgest,p.v.wpplmallong)

combopwpplmal<-fisher(pvalues)


#adult weight and gest length model

cutData <- Data[,c(5,9,10,11,13,38,30,42),drop=FALSE] 
cutData[cutData == -1 ] <-NA
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
wgt.gest<-pglsSEyPagel(Gestation.months.~log10(adult_weight.g.),data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(wgt.gest) 


#grab r squared, lambda, and p values from summary 
r.v.wgt.gest <- summary(wgt.gest)$corBeta
r.v.wgt.gest <- format(r.v.wgt.gest[2,1])
r.v.wgt.gest <-signif(as.numeric(r.v.wgt.gest)^2, digits= 2)
ld.v.wgt.gest<- summary(wgt.gest)$modelStruct$corStruct
ld.v.wgt.gest <- signif(ld.v.wgt.gest[1], digits = 2)
p.v.wgt.gest<-summary(wgt.gest)$tTable
p.v.wgt.gest<-signif(p.v.wgt.gest[2,4], digits = 2)



#plot
wgtgest<-ggplot(cutData, aes(y=Gestation.months., x=log10(adult_weight.g.)))+
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ),)+
  geom_abline() +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Gestation length") +
  xlab("(log10) Adult Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  geom_text_repel(aes(label=ifelse( Gestation.months. > 8,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(title = "Gestation v. Adult Weight",
       subtitle =bquote(p-value:.(p.v.wgt.gest)~R^2:.(r.v.wgt.gest)~Lambda:.(ld.v.wgt.gest))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")


wgtgest


