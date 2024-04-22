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
library(rr2)
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
#adult weight prop

cutData <- Data[,c(5,9,10,11,18,38,42),drop=FALSE] 
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
adult.weight.prop<-pglsSEyPagel(PropMalignant~log10(adult_weight.g.),data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.prop) 



#grab r squared, lambda, and p values from summary 

r.v.adult.weight.prop <- R2(phy = pruned.tree,adult.weight.prop)
r.v.adult.weight.prop <- format(r.v.adult.weight.prop[3])
r.v.adult.weight.prop <-signif(as.numeric(r.v.adult.weight.prop), digits= 2)
ld.v.adult.weight.prop<- summary(adult.weight.prop)$modelStruct$corStruct
ld.v.adult.weight.prop <- signif(ld.v.adult.weight.prop[1], digits = 2)
p.v.adult.weight.prop<-summary(adult.weight.prop)$tTable
p.v.adult.weight.prop<-signif(p.v.adult.weight.prop[2,4], digits = 2)


#plot
ggplot(cutData, aes(y=PropMalignant*100, x=log10(adult_weight.g.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(adult.weight.prop)[1]*100, slope =  coef(adult.weight.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Adult Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. BMR",  
       subtitle =bquote(p-value:.(p.v.adult.weight.prop)~R^2:.(r.v.adult.weight.prop)~Lambda:.(ld.v.adult.weight.prop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.))),max(log10(cutData$adult_weight.g.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=2, y=83.8, label = "62", size = 7)


ggsave(filename='S62wgtpropl.pdf', width=18, height=10, limitsize=FALSE,bg="white")

#gestation models
#gestation prop
cutData <- Data[,c(5,9,10,11,18,30,42),drop=FALSE] 
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
gestation.prop<-pglsSEyPagel(PropMalignant~log10(Gestation.months.),data=cutData,
                            tree=pruned.tree,method="ML",se=SE)

summary(gestation.prop)

#grab r squared, lambda, and p values from summary 

r.v.gestprop <- R2(phy = pruned.tree,gestation.prop)
r.v.gestprop <- format(r.v.gestprop[3])
r.v.gestprop<-signif(as.numeric(r.v.gestprop), digits= 2)
ld.v.gestprop<- summary(gestation.prop)$modelStruct$corStruct
ld.v.gestprop <- signif(ld.v.gestprop[1], digits = 2)
p.v.gestprop<-summary(gestation.prop)$tTable
p.v.gestprop<-signif(p.v.gestprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(Gestation.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(gestation.prop)[1]*100, slope =  coef(gestation.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Gestation (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. BMR",  
       subtitle =bquote(p-value:.(p.v.gestprop)~R^2:.(r.v.gestprop)~Lambda:.(ld.v.gestprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$Gestation.months.))),max(log10(cutData$Gestation.months.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=2, y=83.8, label = "64", size = 7)

ggsave(filename='S64gestprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")


#ggsave(filename='gestpropmal.pdf', width=9.5, height=18, limitsize=FALSE,bg="white")


#litter size models 
#litter size prop
cutData <- Data[,c(5,9,10,11,18,33,42),drop=FALSE] 
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
litter.prop<-pglsSEyPagel(PropMalignant~log10(litter_size),data=cutData,
                         tree=pruned.tree,method="ML",se=SE)
summary(litter.prop)

#grab r squared, lambda, and p values from summary 

r.v.litprop <- R2(phy = pruned.tree,litter.prop)
r.v.litprop <- format(r.v.litprop[3])
r.v.litprop <-signif(as.numeric(r.v.litprop), digits= 2)
ld.v.litprop<- summary(litter.prop)$modelStruct$corStruct
ld.v.litprop <- signif(ld.v.litprop[1], digits = 2)
p.v.litprop<-summary(litter.prop)$tTable
p.v.litprop<-signif(p.v.litprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(litter_size))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(litter.prop)[1]*100, slope =  coef(litter.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Longevity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. BMR",  
       subtitle =bquote(p-value:.(p.v.litprop)~R^2:.(r.v.litprop)~Lambda:.(ld.v.litprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$litter_size))),max(log10(cutData$litter_size))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=2, y=83.8, label = "65", size = 7)


ggsave(filename='S65litprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")


### Longevity model
#longevity prop
cutData <- Data[,c(5,9,10,11,18,40,42),drop=FALSE] 
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
longevity.prop<-pglsSEyPagel(PropMalignant~max_longevity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevity.prop)

#grab r squared, lambda, and p values from summary 

r.v.longprop <- R2(phy = pruned.tree,longevity.prop)
r.v.longprop <- format(r.v.longprop[3])
r.v.longprop <-signif(as.numeric(r.v.longprop), digits= 2)
ld.v.longprop<- summary(longevity.prop)$modelStruct$corStruct
ld.v.longprop <- signif(ld.v.longprop[1], digits = 2)
p.v.longprop<-summary(longevity.prop)$tTable
p.v.longprop<-signif(p.v.longprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(max_longevity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(longevity.prop)[1]*100, slope =  coef(longevity.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Longevity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. BMR",  
       subtitle =bquote(p-value:.(p.v.longprop)~R^2:.(r.v.longprop)~Lambda:.(ld.v.longprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$max_longevity.months.))),max(log10(cutData$max_longevity.months.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=2, y=83.8, label = "66", size = 7)


ggsave(filename='S66longprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")

##BMR models
#bmr prop
cutData <- Data[,c(5,9,10,11,18,41,42),drop=FALSE] 
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
BMR.prop<-pglsSEyPagel(PropMalignant~log10(metabolic_rate),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)
summary(BMR.prop)

#grab r squared, lambda, and p values from summary 

r.v.bmrprop <- R2(phy = pruned.tree,BMR.prop)
r.v.bmrprop <- format(r.v.bmrprop[3])
r.v.bmrprop <-signif(as.numeric(r.v.bmrprop), digits= 2)
ld.v.bmrprop<- summary(BMR.prop)$modelStruct$corStruct
ld.v.bmrprop <- signif(ld.v.bmrprop[1], digits = 2)
p.v.bmrprop<-summary(BMR.prop)$tTable
p.v.bmrprop<-signif(p.v.bmrprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(metabolic_rate))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(BMR.prop)[1]*100, slope =  coef(BMR.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) basal metabolic rate") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. BMR",  
       subtitle =bquote(p-value:.(p.v.bmrprop)~R^2:.(r.v.bmrprop)~Lambda:.(ld.v.bmrprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$metabolic_rate))),max(log10(cutData$metabolic_rate))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=2, y=83.8, label = "67", size = 7)

ggsave(filename='S67bmrprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")


#wxl models #weight and longevity
#wxl prop
cutData <- Data[,c(5,9,10,11,18,40,38,42),drop=FALSE] 
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
wxl.prop<-pglsSEyPagel(PropMalignant~log10(max_longevity.months.*adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wxl.prop)

#grab r squared, lambda, and p values from summary 

r.v.wxprop <- R2(phy = pruned.tree,wxl.prop)
r.v.wxprop <- format(r.v.wxprop[3])
r.v.wxprop <-signif(as.numeric(r.v.wxprop), digits= 2)
ld.v.wxprop<- summary(wxl.prop)$modelStruct$corStruct
ld.v.wxprop <- signif(ld.v.wxprop[1], digits = 2)
p.v.wxprop<-summary(wxl.prop)$tTable
p.v.wxprop<-signif(p.v.wxprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(adult_weight.g.*max_longevity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(wxl.prop)[1]*100, slope =  coef(wxl.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Adult Weight(g)*Longevity(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Max Longevity*Weight",  
       subtitle =bquote(p-value:.(p.v.wxprop)~R^2:.(r.v.wxprop)~Lambda:.(ld.v.wxprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.*cutData$max_longevity.months.))),max(log10(cutData$adult_weight.g.*cutData$max_longevity.months.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=2, y=83.8, label = "68", size = 7)

ggsave(filename='S68wgtlongprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")

#litters per year prop
cutData <- Data[,c(5,9,10,11,18,34,42),drop=FALSE] 
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
lityear.prop<-pglsSEyPagel(PropMalignant~log10(litters_year),data=cutData,
                          tree=pruned.tree,method="ML",se=SE)
summary(lityear.prop)

#grab r squared, lambda, and p values from summary 

r.v.lyearprop <- R2(phy = pruned.tree,lityear.prop)
r.v.lyearprop <- format(r.v.lyearprop[3])
r.v.lyearprop<-signif(as.numeric(r.v.lyearprop), digits= 2)
ld.v.lyearprop<- summary(lityear.prop)$modelStruct$corStruct
ld.v.lyearprop <- signif(ld.v.lyearprop[1], digits = 2)
p.v.lyearprop<-summary(lityear.prop)$tTable
p.v.lyearprop<-signif(p.v.lyearprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(litters_year))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(GrowthR.prop)[1]*100, slope =  coef(GrowthR.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Litters Per YEar") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse( PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Litters per Year",  
       subtitle =bquote(p-value:.(p.v.lyearprop)~R^2:.(r.v.lyearprop)~Lambda:.(ld.v.lyearprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$litters_year)),log10(max(cutData$litters_year))),
                  ylim = c(0,100),clip = "off")+
  annotate("text",x=-3.36, y=83.8, label = "73", size = 7)

ggsave(filename='S69lityearprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")


####Female Maturity prop
cutData <- Data[,c(5,9,10,11,18,28,42),drop=FALSE] 
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
Fmaturity.prop<-pglsSEyPagel(PropMalignant~female_maturity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(Fmaturity.prop)

#grab r squared, lambda, and p values from summary 

r.v.fmaturityprop <- R2(phy = pruned.tree,Fmaturity.prop)
r.v.fmaturityprop <- format(r.v.fmaturityprop[3])
r.v.fmaturityprop<-signif(as.numeric(r.v.fmaturityprop), digits= 2)
ld.v.fmaturityprop<- summary(Fmaturity.prop)$modelStruct$corStruct
ld.v.fmaturityprop <- signif(ld.v.fmaturityprop[1], digits = 2)
p.v.fmaturityprop<-summary(Fmaturity.prop)$tTable
p.v.fmaturityprop<-signif(p.v.fmaturityprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(female_maturity.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(Fmaturity.prop)[1]*100, slope =  coef(Fmaturity.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("Female Maturity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse( PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Female Maturity",  
       subtitle =bquote(p-value:.(p.v.fmaturityprop)~R^2:.(r.v.fmaturityprop)~Lambda:.(ld.v.fmaturityprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$female_maturity.months.)),log10(max(cutData$female_maturity.months.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=-.09, y=83.9, label = "70", size = 7)

ggsave(filename='S70femmatprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")


####Male Maturity prop
cutData <- Data[,c(5,9,10,11,18,29,42),drop=FALSE] 
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
Mmaturity.prop<-pglsSEyPagel(PropMalignant~male_maturity.months.,data=cutData,
                            tree=pruned.tree,method="ML",se=SE)
summary(Mmaturity.prop)

#grab r squared, lambda, and p values from summary 

r.v.Mmaturityprop <- R2(phy = pruned.tree,Mmaturity.prop)
r.v.Mmaturityprop <- format(r.v.Mmaturityprop[3])
r.v.Mmaturityprop<-signif(as.numeric(r.v.Mmaturityprop), digits= 2)
ld.v.Mmaturityprop<- summary(Mmaturity.prop)$modelStruct$corStruct
ld.v.Mmaturityprop <- signif(ld.v.Mmaturityprop[1], digits = 2)
p.v.Mmaturityprop<-summary(Mmaturity.prop)$tTable
p.v.Mmaturityprop<-signif(p.v.Mmaturityprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(male_maturity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(Mmaturity.prop)[1]*100, slope =  coef(Mmaturity.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("Male Maturity (Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Male Maturity",  
       subtitle =bquote(p-value:.(p.v.Mmaturityprop)~R^2:.(r.v.Mmaturityprop)~Lambda:.(ld.v.Mmaturityprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$male_maturity.months.)),log10(max(cutData$male_maturity.months.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=-0.74, y=83.8, label = "71", size = 7)

ggsave(filename='S71malematprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")

#Weaning Weight prop
cutData <- Data[,c(5,9,10,11,18,37,42),drop=FALSE] 
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
weanw.prop<-pglsSEyPagel(PropMalignant~log10(weaning_weight.g.),data=cutData,
                        tree=pruned.tree,method="ML",se=SE)
summary(weanw.prop)

#grab r squared, lambda, and p values from summary 

r.v.weanwprop <- R2(phy = pruned.tree,weanw.prop)
r.v.weanwprop <- format(r.v.weanwprop[3])
r.v.weanwprop<-signif(as.numeric(r.v.weanwprop), digits= 2)
ld.v.weanwprop<- summary(weanw.prop)$modelStruct$corStruct
ld.v.weanwprop <- signif(ld.v.weanwprop[1], digits = 2)
p.v.weanwprop<-summary(weanw.prop)$tTable
p.v.weanwprop<-signif(p.v.weanwprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(weaning_weight.g.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(weanw.prop)[1]*100, slope =  coef(weanw.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Weaning Weight (g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse( PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Weaning Weight",  
       subtitle =bquote(p-value:.(p.v.weanwprop)~R^2:.(r.v.weanwprop)~Lambda:.(ld.v.weanwprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$weaning_weight.g.)),log10(max(cutData$weaning_weight.g.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text",  x=-.64, y=83.8, label = "72", size = 7)

ggsave(filename='S72weanprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")

#Growth Rate prop
cutData <- Data[,c(5,9,10,11,18,39,42),drop=FALSE] 
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
GrowthR.prop<-pglsSEyPagel(PropMalignant~log10(growth_rate.1.days.),data=cutData,
                          tree=pruned.tree,method="ML",se=SE)
summary(GrowthR.prop)

#grab r squared, lambda, and p values from summary 

r.v.GrowthRprop <- R2(phy = pruned.tree,GrowthR.prop)
r.v.GrowthRprop <- format(r.v.GrowthRprop[3])
r.v.GrowthRprop<-signif(as.numeric(r.v.GrowthRprop), digits= 2)
ld.v.GrowthRprop<- summary(GrowthR.prop)$modelStruct$corStruct
ld.v.GrowthRprop <- signif(ld.v.GrowthRprop[1], digits = 2)
p.v.GrowthRprop<-summary(GrowthR.prop)$tTable
p.v.GrowthRprop<-signif(p.v.GrowthRprop[2,4], digits = 3)

ggplot(cutData, aes(y=PropMalignant*100, x=log10(growth_rate.1.days.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(GrowthR.prop)[1]*100, slope =  coef(GrowthR.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Growth Rate") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse( PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Growth Rate",  
       subtitle =bquote(p-value:.(p.v.GrowthRprop)~R^2:.(r.v.GrowthRprop)~Lambda:.(ld.v.GrowthRprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c(log10(min(cutData$growth_rate.1.days.)),log10(max(cutData$growth_rate.1.days.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text",x=-3.36, y=83.8, label = "73", size = 7)

ggsave(filename='S73growprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")

#w+g

cutData <- Data[,c(5,9,10,11,18,38,30,42),drop=FALSE] 
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
wpl.prop<-pglsSEyPagel(PropMalignant~log10(Gestation.months.)+log10(adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

coef(wpl.prop)

summary(wpl.prop)

#grab r squared, lambda, and p values from summary 

r.v.wpprop <- R2(phy = pruned.tree,wpl.prop)
r.v.wpprop <- format(r.v.wpprop[3])
r.v.wpprop <-signif(as.numeric(r.v.wpprop), digits= 2)
ld.v.wpprop<- summary(wpl.prop)$modelStruct$corStruct
ld.v.wpprop <- signif(ld.v.wpprop[1], digits = 2)
p.v.wpprop<-summary(wpl.prop)$tTable
p.v.wppropgest<-signif(p.v.wpprop[2,4], digits = 3)
p.v.wppropwgt<-signif(p.v.wpprop[3,4], digits = 2)
p.v.wpgpropgest<-signif(p.v.wpprop[2,4], digits = 3)
p.v.wpgpropgest<-signif(p.v.wpprop[3,4], digits = 2)
c.wpplpropgest<-coef(wpl.prop)[2]
c.wpplproplong<-coef(wpl.prop)[3]


# Assuming you have a vector of p-values from multiple tests
p.values <- c(p.v.wppropgest,p.v.wppropwgt)

# Use the p.adjust function with the BH method to adjust the p-values
adjusted.p.values <- p.adjust(p.values, method = "BH")

# Determine which coefficients are significant at the 10% FDR level
significant_effects <- adjusted.p.values < 0.1

# Print which coefficients are significant
print(significant_effects)



ggplot(cutData, aes(y=PropMalignant*100, x=log10(adult_weight.g.)+log10(Gestation.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(wpl.prop)[1]*100, slope =  coef(wpl.prop)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Prop Malignant (%)") +
  xlab("(log10) Adult Weight(g)+Gestation(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  #geom_text_repel(aes(label=ifelse(PropMalignant > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Prop Malignant vs. Max Longevity+Gestation",  
       subtitle =bquote(p-value:.(p.v.wpprop)~R^2:.(r.v.wpprop)~Lambda:.(ld.v.wpprop))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),
                  ylim = c(0,100),clip = "off")+
  annotate("text", x=.18, y=83.8, label = "74", size = 7)

ggsave(filename='S74wgtgestprop.pdf', width=18, height=10, limitsize=FALSE,bg="white")



#g+l

cutData <- Data[,c(5,9,10,11,18,40,30,42),drop=FALSE] 
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
wppl.prop<-pglsSEyPagel(PropMalignant~log10(Gestation.months.)+log10(max_longevity.months.),data=cutData,
                       tree=pruned.tree,method="ML",se=SE)

summary(wppl.prop)
coef(wppl.prop)

#grab r squared, lambda, and p values from summary 

r.v.wpplprop <- R2(phy = pruned.tree,wppl.prop)
r.v.wpplprop <- format(r.v.wpplprop[3])
r.v.wpplprop <-signif(as.numeric(r.v.wpplprop), digits= 2)
ld.v.wpplprop<- summary(wppl.prop)$modelStruct$corStruct
ld.v.wpplprop <- signif(ld.v.wpplprop[1], digits = 2)
p.v.wpplprop<-summary(wppl.prop)$tTable
p.v.wpplpropgest<-signif(p.v.wpplprop[2,4], digits = 3)
p.v.wpplproplong<-signif(p.v.wpplprop[3,4], digits = 2)
c.wpplpropgest<-coef(wppl.prop)[2]
c.wpplproplong<-coef(wppl.prop)[3]


