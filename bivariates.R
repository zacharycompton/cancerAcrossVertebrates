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
Data <- read.csv("min20516.csv")
View(Data)


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

wxl.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(max_longevity.months.*adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wxl.neo)

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

ggsave(filename='wgtlongneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

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

wxl.mal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity.months.*adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wxl.mal)

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

ggsave(filename='wgtlongmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

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

wpl.neo<-pglsSEyPagel(NeoplasiaPrevalence~Gestation.months.+log10(adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.neo)

r.v.wpneo <- summary(wpl.neo)$corBeta
r.v.wpneo <- format(r.v.wpneo[2,1])
r.v.wpneo <-signif(as.numeric(r.v.wpneo)^2, digits= 2)
ld.v.wpneo<- summary(wpl.neo)$modelStruct$corStruct
ld.v.wpneo <- signif(ld.v.wpneo[1], digits = 2)
p.v.wpneo<-summary(wpl.neo)$tTable
p.v.wpneo<-signif(p.v.wpneo[2,4], digits = 3)

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

ggsave(filename='gplongneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")


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

wpl.mal<-pglsSEyPagel(MalignancyPrevalence~Gestation.months.+log10(adult_weight.g.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.mal)

r.v.wpmal <- summary(wpl.mal)$corBeta
r.v.wpmal <- format(r.v.wpmal[2,1])
r.v.wpmal <-signif(as.numeric(r.v.wpmal)^2, digits= 2)
ld.v.wpmal<- summary(wpl.mal)$modelStruct$corStruct
ld.v.wpmal <- signif(ld.v.wpmal[1], digits = 2)
p.v.wpmal<-summary(wpl.mal)$tTable
p.v.wpmal<-signif(p.v.wpmal[2,4], digits = 3)

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


ggsave(filename='gplongmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")





#gestm + litter size

cutData <- Data[,c(5,9,10,11,13,38,33,42),drop=FALSE] 
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

lpg.neo<-pglsSEyPagel(NeoplasiaPrevalence~Gestation.months.+log10(litter_size),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(lpg.neo)

r.v.lpgneo <- summary(lpg.neo)$corBeta
r.v.lpgneo <- format(r.v.lpgneo[2,1])
r.v.lpgneo <-signif(as.numeric(r.v.lpgneo)^2, digits= 2)
ld.v.lpgneo<- summary(lpg.neo)$modelStruct$corStruct
ld.v.lpgneo <- signif(ld.v.lpgneo[1], digits = 2)
p.v.lpgneo<-summary(lpg.neo)$tTable
p.v.lpgneo<-signif(p.v.lpgneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(litter_size)+log10(Gestation.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(lpg.neo)[1]*100, slope =  coef(lpg.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_colpgot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Litter Size+Gestation(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence +. Litter Size+Gestation",  
       subtitle =bquote(p-value:.(p.v.lpgneo)~R^2:.(r.v.lpgneo)~Lambda:.(ld.v.lpgneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$litter_size)+log10(cutData$Gestation.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=.13, y=83.8, label = "19", size = 7)

ggsave(filename='gplitneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")


#lit+g mal

cutData <- Data[,c(5,9,10,11,17,38,33,42),drop=FALSE] 
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

wpl.mal<-pglsSEyPagel(MalignancyPrevalence~Gestation.months.+log10(litter_size),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.mal)

r.v.wpmal <- summary(wpl.mal)$corBeta
r.v.wpmal <- format(r.v.wpmal[2,1])
r.v.wpmal <-signif(as.numeric(r.v.wpmal)^2, digits= 2)
ld.v.wpmal<- summary(wpl.mal)$modelStruct$corStruct
ld.v.wpmal <- signif(ld.v.wpmal[1], digits = 2)
p.v.wpmal<-summary(wpl.mal)$tTable
p.v.wpmal<-signif(p.v.wpmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(Gestation.months.)+log10(litter_sie))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c((min(log10(cutData$litter_size)+log10(cutData$Gestation.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$Gestation.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(wpl.mal)[1]*100, slope =  coef(wpl.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Litter Size+Gestation(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Litter Sizey+Gestation",  
       subtitle =bquote(p-value:.(p.v.wpmal)~R^2:.(r.v.wpmal)~Lambda:.(ld.v.wpmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")+
  annotate("text", x=.13, y=50.3, label = "20", size = 7)


ggsave(filename='gplitmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")



#wpl neo
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

wpl.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(adult_weight.g.)+log10(max_longevity.months.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.neo)

r.v.wplneo <- summary(wpl.neo)$corBeta
r.v.wplneo <- format(r.v.wplneo[2,1])
r.v.wplneo <-signif(as.numeric(r.v.wplneo)^2, digits= 2)
ld.v.wplneo<- summary(wpl.neo)$modelStruct$corStruct
ld.v.wplneo <- signif(ld.v.wplneo[1], digits = 2)
p.v.wplneo<-summary(wpl.neo)$tTable
p.v.wplneo<-signif(p.v.wplneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(adult_weight.g.)+log10(max_longevity.months.))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(wpl.neo)[1]*100, slope =  coef(wpl.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Adult Weight(g)+Longevity(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Max Longevity+Weight",  
       subtitle =bquote(p-value:.(p.v.wplneo)~R^2:.(r.v.wplneo)~Lambda:.(ld.v.wplneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.)+log10(cutData$max_longevity.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$max_longevity.months.))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=2, y=83.8, label = "8", size = 7)

ggsave(filename='wplongneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#wpl mal
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

wpl.mal<-pglsSEyPagel(MalignancyPrevalence~log10(adult_weight.g.)+log10(max_longevity.months.),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.mal)

r.v.wplmal <- summary(wpl.mal)$corBeta
r.v.wplmal <- format(r.v.wplmal[2,1])
r.v.wplmal <-signif(as.numeric(r.v.wplmal)^2, digits= 2)
ld.v.wplmal<- summary(wpl.mal)$modelStruct$corStruct
ld.v.wplmal <- signif(ld.v.wplmal[1], digits = 2)
p.v.wplmal<-summary(wpl.mal)$tTable
p.v.wplmal<-signif(p.v.wplmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(adult_weight.g.)+log10(max_longevity.months.))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c((min(log10(cutData$adult_weight.g.)+log10(cutData$max_longevity.months.))),max(log10(cutData$adult_weight.g.)+log10(cutData$max_longevity.months.))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(wpl.mal)[1]*100, slope =  coef(wpl.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Adult Weight(g)+Longevity(Mo)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Max Longevity+Weight",  
       subtitle =bquote(p-value:.(p.v.wplmal)~R^2:.(r.v.wplmal)~Lambda:.(ld.v.wplmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  annotate("text", x=2, y=50.3, label = "9", size = 7)

ggsave(filename='wplongmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#wplit neo
cutData <- Data[,c(5,9,10,11,13,33,38,42),drop=FALSE] 
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

wpl.neo<-pglsSEyPagel(NeoplasiaPrevalence~log10(litter_size)+log10(adult_weight.g),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.neo)

r.v.wplitneo <- summary(wpl.neo)$corBeta
r.v.wplitneo <- format(r.v.wplitneo[2,1])
r.v.wplitneo <-signif(as.numeric(r.v.wplitneo)^2, digits= 2)
ld.v.wplitneo<- summary(wpl.neo)$modelStruct$corStruct
ld.v.wplitneo <- signif(ld.v.wplitneo[1], digits = 2)
p.v.wplitneo<-summary(wpl.neo)$tTable
p.v.wplitneo<-signif(p.v.wplitneo[2,4], digits = 3)

ggplot(cutData, aes(y=NeoplasiaPrevalence*100, x=log10(litter_size)+log10(adult_weight.g))) + 
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(wpl.neo)[1]*100, slope =  coef(wpl.neo)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Neoplasia Prevalence (%)") +
  xlab("(log10) Adult Weight(g)+Litter Size") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(NeoplasiaPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Neoplasia Prevalence vs. Litter Size+Weight",  
       subtitle =bquote(p-value:.(p.v.wplitneo)~R^2:.(r.v.wplitneo)~Lambda:.(ld.v.wplitneo))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  guides(size=guide_legend())+
  coord_cartesian(xlim = c((min(log10(cutData$litter_size)+log10(cutData$adult_weight.g))),max(log10(cutData$litter_size)+log10(cutData$adult_weight.g))),
                  ylim = c(0,75),clip = "off")+
  annotate("text", x=2, y=83.8, label = "8", size = 7)

ggsave(filename='wplitneo.pdf', width=13, height=10, limitsize=FALSE,bg="white")

#wplit mal
cutData <- Data[,c(5,9,10,11,17,33,38,42),drop=FALSE] 
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

wpl.mal<-pglsSEyPagel(MalignancyPrevalence~log10(litter_size)+log10(adult_weight.g),data=cutData,
                      tree=pruned.tree,method="ML",se=SE)

summary(wpl.mal)

r.v.wplitmal <- summary(wpl.mal)$corBeta
r.v.wplitmal <- format(r.v.wplitmal[2,1])
r.v.wplitmal <-signif(as.numeric(r.v.wplitmal)^2, digits= 2)
ld.v.wplitmal<- summary(wpl.mal)$modelStruct$corStruct
ld.v.wplitmal <- signif(ld.v.wplitmal[1], digits = 2)
p.v.wplitmal<-summary(wpl.mal)$tTable
p.v.wplitmal<-signif(p.v.wplitmal[2,4], digits = 3)

ggplot(cutData, aes(y=MalignancyPrevalence*100, x=log10(litter_size)+log10(adult_weight.g))) +
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff"))+
  scale_y_continuous(
    breaks = c(0, 25,45),
    labels = c(0, 25,45))+
  coord_cartesian(xlim = c((min(log10(cutData$litter_size)+log10(cutData$adult_weight.g))),max(log10(cutData$litter_size)+log10(cutData$adult_weight.g))),
                  ylim = c(0,45),clip = "off")+
  geom_abline(intercept = coef(wpl.mal)[1]*100, slope =  coef(wpl.mal)[2]*100,
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("(log10) Litter Size+Adult Weight(g)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  geom_text_repel(aes(label=ifelse(MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  labs(title = "Malignancy Prevalence vs. Litter Size+Weight",  
       subtitle =bquote(p-value:.(p.v.wplitmal)~R^2:.(r.v.wplitmal)~Lambda:.(ld.v.wplitmal))) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+   labs(colour="Clade", size="Total Necropsies")+
  annotate("text", x=2, y=50.3, label = "9", size = 7)

ggsave(filename='wplitmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")


