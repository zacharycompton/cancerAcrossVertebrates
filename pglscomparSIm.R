library(ape)
library(nlme)
library(caper)
library(dplyr)

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
Data <- read.csv("min20-2022.05.16.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
View(Data)


#adult weight models
#adult weight neo

cutData <- Data[,c(5,9,10,11,12,13,38,42),drop=FALSE] 
cutData[cutData < 0] <-NA
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

generate_strictly_positive_data <- function(n, mean, sd) {
  data <- rnorm(n, mean, sd) # Generate data
  while(any(data <= 0)) { # Check for non-positive values
    data[data <= 0] <- rnorm(sum(data <= 0), mean, sd) # Re-sample non-positives
  }
  return(data)
}

i=0
resultsAWSim<-data.frame(pvalPGLS = NULL, slopePGLS = NULL, pvalGEE = NULL, slopeGEE= NULL)
for (i in i:100){
# NeoplasiaPrevSim<-rlnorm(nrow(cutData), mean(cutData$NeoplasiaPrevalence), sd(cutData$NeoplasiaPrevalence))
AdultWSim<-generate_strictly_positive_data(nrow(cutData), mean(cutData$adult_weight.g.), sd(cutData$adult_weight.g.))

data_simulated <- data.frame(TotalRecords= cutData$RecordsWithDenominators,NeoplasiaPrev = cutData$NeoplasiaPrevalence, AdultWSim = AdultWSim, SE = cutData$SE_simple, NeoplasiaOccurences = cutData$NeoplasiaWithDenominators, NonOccurences = cutData$RecordsWithDenominators-cutData$NeoplasiaWithDenominators, row.names = row.names(cutData))

NeoplasiaOccurences<-data_simulated$NeoplasiaOccurences
NonOccurences<-data_simulated$NonOccurences


#pgls model
adult.weight.neo<-pglsSEyPagel(NeoplasiaPrev~AdultWSim,data=data_simulated,tree=pruned.tree,se=SE,method = "ML")

summary(adult.weight.neo) 


#compar 
compar<-invisible(compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~AdultWSim, phy = pruned.tree, family = "binomial"))


#grab r squared, lambda, and p values from summary 
p.v.adult.weight.neo<-summary(adult.weight.neo)$tTable
p.v.adult.weight.neo<-signif(p.v.adult.weight.neo[2,4], digits = 2)


geePval<-foo(compar)[2,4]
geeSlope<-foo(compar)[2,1]

new_row <- data.frame(pvalPGLS = p.v.adult.weight.neo, slopePGLS = adult.weight.neo$coefficients[2], pvalGEE = geePval, slopeGEE=geeSlope )

# Add the new row at the end of the dataframe
resultsAWSim <- rbind(resultsAWSim, new_row)


}

resultsLongSim<-data.frame(pvalPGLS = NULL, slopePGLS = NULL, pvalGEE = NULL, slopeGEE= NULL)
i=0



cutData <- Data[,c(5,9,10,11,12,13,40,42),drop=FALSE] 
cutData[cutData < 0] <-NA
cutData <- na.omit(cutData)
tree <- read.tree("min20Fixed516.nwk")


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

for (i in i:100){
### Longevity model
#longevity neo

LongSim<-generate_strictly_positive_data(nrow(cutData), mean(cutData$max_longevity.months.), sd(cutData$max_longevity.months.))


data_simulated <- data.frame(TotalRecords= cutData$RecordsWithDenominators,NeoplasiaPrev = cutData$NeoplasiaPrevalence, LongSim = LongSim, SE = cutData$SE_simple, NeoplasiaOccurences = cutData$NeoplasiaWithDenominators, NonOccurences = cutData$RecordsWithDenominators-cutData$NeoplasiaWithDenominators, row.names = row.names(cutData))

NeoplasiaOccurences<-data_simulated$NeoplasiaOccurences
NonOccurences<-data_simulated$NonOccurences


longevity.neo<-pglsSEyPagel(NeoplasiaPrev~LongSim,data=data_simulated,
                            tree=pruned.tree,method="ML",se=SE)
summary(longevity.neo)

compar<-invisible(compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~LongSim, phy = pruned.tree, family = "binomial"))


geePval<-foo(compar)[2,4]
geeSlope<-foo(compar)[2,1]



#grab r squared, p value, and lambda from summary 

r.v.longneo <- summary(longevity.neo)$corBeta
r.v.longneo <- format(r.v.longneo[2,1])
r.v.longneo <-signif(as.numeric(r.v.longneo)^2, digits= 2)
ld.v.longneo<- summary(longevity.neo)$modelStruct$corStruct
ld.v.longneo <- signif(ld.v.longneo[1], digits = 2)
p.v.longneo<-summary(longevity.neo)$tTable
p.v.longneo<-signif(p.v.longneo[2,4], digits = 3)



geePval<-foo(compar)[2,4]
geeSlope<-foo(compar)[2,1]

new_row <- data.frame(pvalPGLS = p.v.longneo, slopePGLS = longevity.neo$coefficients[2], pvalGEE = geePval, slopeGEE=geeSlope )

# Add the new row at the end of the dataframe
resultsLongSim <- rbind(resultsLongSim, new_row)



}

write.csv(resultsLongSim, "./longSimResults.csv")
write.csv(resultsAWSim, "./adultWeightResults.csv")




