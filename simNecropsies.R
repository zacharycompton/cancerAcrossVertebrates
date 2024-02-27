library(MASS)
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
library(lme4)

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



data<-read.csv("min1LH.csv")
# Fit a Negative Binomial model to the 'Necropsies' column
necropsy_counts <- data$RecordsWithDenominators
fit <- glm.nb(necropsy_counts ~ 1) # fitting to an intercept only model

# Summary of the fitted model
summary(fit)

# Simulate new data based on the fitted model
# The number to simulate, `n`, should be set to the number of species you want to simulate data for
data$Prob_Neoplasia <- ifelse(data$RecordsWithDenominators > 0, 
                                       data$NeoplasiaWithDenominators / data$RecordsWithDenominators, 
                                       0)

# Simulate neoplasias for each species
simulated_neoplasias <- rbinom(n, size = simulated_necropsies, prob = data$Prob_Neoplasia)

# Ensure simulated neoplasias are not greater than simulated necropsies
simulated_neoplasias <- pmin(simulated_neoplasias, simulated_necropsies)


# Add the simulated values to the data frame
data$SimulatedNecropsies <- simulated_necropsies
data$SimulatedNeoplasias <- simulated_neoplasias
data$SE_SimpleSim <- ifelse(data$RecordsWithDenominators > 0, 
                                  1 / data$RecordsWithDenominators, 
                                  NA)
data$NeoplasiaPrevalenceSim<-data$SimulatedNeoplasias/data$SimulatedNecropsies 


species_counts <- data %>%
  group_by(SimulatedNecropsies) %>%
  summarise(SpeciesCount = n_distinct(Species))


ggplot(species_counts, aes(x = SimulatedNecropsies, y = SpeciesCount)) +
  geom_col() +  # geom_col() is used for bar plots, which this effectively is
  theme_cowplot(12) +
  labs(x = "Number of Records", y = "Number of Species", title = "Distribution of Species by Number of Records")

neoplasia_count <- data %>%
  group_by(SimulatedNeoplasias) %>%
  summarise(SpeciesCount = n_distinct(Species))


ggplot(neoplasia_count, aes(x = SimulatedNeoplasias, y = SpeciesCount)) +
  geom_col() +  # geom_col() is used for bar plots, which this effectively is
  theme_cowplot(12) +
  labs(x = "Number of Neoplasias", y = "Number of Species", title = "Distribution of Species by Number of Records")



real_neoplasia_count <- data %>%
  group_by(NeoplasiaWithDenominators) %>%
  summarise(SpeciesCount = n_distinct(Species))


ggplot(real_neoplasia_count, aes(x = NeoplasiaWithDenominators, y = SpeciesCount)) +
  geom_col() +  # geom_col() is used for bar plots, which this effectively is
  theme_cowplot(12) +
  labs(x = "Number of Neoplasias Real", y = "Number of Species", title = "Distribution of Species by Number of Records")



#adult weight models
#adult weight neo

cutData <- data[,c(5,9,11,13,38,42,43,44,45,46),drop=FALSE] 
condition <- cutData$adult_weight == -1 & !is.na(cutData$adult_weight)
cutData[condition, ] <- NA

cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")

cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_SimpleSim,cutData$Species)[rownames(cutData)]
view(cutData)


#pgls model
adult.weight.neo<-pglsSEyPagel(NeoplasiaPrevalenceSim~log10(adult_weight),data=cutData,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.neo) 

#grab r squared, lambda, and p values from summary 

r.v.adult.weight.neo <- summary(adult.weight.neo)$corBeta
r.v.adult.weight.neo <- format(r.v.adult.weight.neo[2,1])
r.v.adult.weight.neo <-signif(as.numeric(r.v.adult.weight.neo)^2, digits= 2)
ld.v.adult.weight.neo<- summary(adult.weight.neo)$modelStruct$corStruct
ld.v.adult.weight.neo <- signif(ld.v.adult.weight.neo[1], digits = 2)
p.v.adult.weight.neo<-summary(adult.weight.neo)$tTable
p.v.adult.weight.neo<-signif(p.v.adult.weight.neo[2,4], digits = 2)




adultWeight<-log10(cutData$adult_weight)
NeoplasiaOccurences<-cutData$SimulatedNeoplasias
NonOccurences<-(cutData$SimulatedNecropsies-cutData$SimulatedNeoplasias)

compar<-compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~adultWeight, phy = pruned.tree, family = "binomial")
compar
