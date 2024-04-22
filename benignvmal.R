#open libraries
library(patchwork)
library(tidyverse)
library(cowplot)
library(ggrepel)


#read csv
cleanPath <- read.csv("cleanPath.min20.062822.csv")
min20Data <- read.csv("min20-2022.05.16.csv")
view(min20Data)

#linear model mav vs benign
malben<-lm(MalignancyPrevalence~BenignPrevalence,data = min20Data) 
summary(malben)
summary(malben)$r.squared
summary(malben)$coefficients[2,4]

#ks test for benign lifespan vs malignant lifespan
benign<-filter(cleanPath, is.element(Malignant, c(0)))
malignant<-filter(cleanPath, is.element(Malignant, c(1)))
ks.test(benign$proportion_lifespan,malignant$proportion_lifespan)
cleanPath<- filter(cleanPath,Malignant >= 0)


#grab r squared
r2<-signif(as.numeric(summary(malben)$r.squared)^2, digits= 2)
pval<-signif(as.numeric(summary(malben)$coefficients[2,4])^2, digits= 2)



#plot malignancy prev vs benign prev
ggplot(min20Data, aes(y=MalignancyPrevalence*100, x= BenignPrevalence*100))+
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ),)+
  scale_y_continuous(
    limits = c(0,50),
    breaks = c(0, 25,50),
    labels = c(0, 25,50))+
  scale_x_continuous(
    limits = c(0,30),
    breaks = c(0, 15,30),
    labels = c(0, 15,30)
  )+
  geom_abline(intercept = coef(malben)[1], slope =  coef(malben)[2],
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("Benign Prevalence (%)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(title = "Malignancy Prevalence v. Benign Prevalence",
       subtitle =bquote(p-value:.(pval)~R^2:.(r2)))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")



#save 
ggsave(filename='benignvmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")



#plot prop lifespan malignant factor
bvm<- ggplot(cleanPath, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Malignant)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Slow Life History",
          subtitle="All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c( "No Tumor Found","Benign", "Malignant"))

bvm






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
#adult weight neo

cutData <- Data[,c(5,9,10,11,17,20,42),drop=FALSE] 
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
benmal<-pglsSEyPagel(MalignancyPrevalence~BenignPrevalence,data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(benmal) 



#grab r squared, lambda, and p values from summary 

r.v.benmal <- R2(phy = pruned.tree,benmal)
r.v.benmal <- format(r.v.benmal[3])
r.v.benmal <-signif(as.numeric(r.v.benmal), digits= 2)
ld.v.benmal<- summary(benmal)$modelStruct$corStruct
ld.v.benmal <- signif(ld.v.benmal[1], digits = 2)
p.v.benmal<-summary(benmal)$tTable
p.v.benmal<-signif(p.v.benmal[2,4], digits = 2)




#plot
bem<-ggplot(cutData, aes(y=MalignancyPrevalence*100, x=BenignPrevalence))+
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ),)+
  scale_y_continuous(
    limits = c(0,75),
    breaks = c(0, 25,50,75),
    labels = c(0, 25,50,75))+
  geom_abline(intercept = coef(benmal)[1]*100, slope =  coef(benmal)[2]*100,
              color = 'grey',size = 1.2) +
  labs(title = "Malignancy Prevalence vs. Benign Prevalence",  
       subtitle =bquote(p-value:.(p.v.benmal)~R^2:.(r.v.benmal)~Lambda:.(ld.v.benmal))) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("Benign Prevalence (%)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")

bem


ggsave(filename='S57mutation.mal.png', width=9.5, height=7, limitsize=FALSE,bg="white")
