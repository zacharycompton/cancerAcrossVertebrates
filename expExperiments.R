library(MASS)
library(nlme)
#library(rms)
library(phytools)
library(geiger)
library(caper)
library(tidyverse)
library(parallel)


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


dataReport <- data.frame(n = numeric(),min = numeric(), model = character(), pval = numeric(), simslope = numeric(),slope = numeric())
Data <- read.csv("min1LH.csv") # Adjust filename/path as necessary

# Define parameter settings
min_necropsies_values <- c(5,20)
#min_necropsies_values <- 20
num_species_values <- c(50, 150,300)
#num_species_values <- 50
num_simulations <- 100
tree <- read.tree("min1.nwk")


slope_values <- c(.021, 0, -.002)  # Array of slope values


no_cores <- detectCores() - 1



run_simulation <- function(params) {
  compar_row <- data.frame(n = NA,min = NA, model = "Compar", pval = NA,simslope= NA,  slope = NA)
  # Extract the slope directly from params as it's directly passed
  slope <- params
  
  # Initialize an empty list to store results of each simulation
  all_results <- list()
  
  for(min_necropsies in min_necropsies_values) {
    for(num_species in num_species_values) {
      for(sim in 1:num_simulations) {
      
# Read Data
      # Sample 'num_species' species from the dataset
      mamData <- Data[Data$RecordsWithDenominators >= min_necropsies & Data$Class == "Mammalia", ]
      mamData <- mamData[sample(nrow(mamData), min(num_species, nrow(mamData))), ]
      n <- nrow(mamData)
      maxNeoPrev<-max(na.omit(mamData$NeoplasiaPrevalence))
      maxNeo<-max(mamData$NeoplasiaWithDenominators)
      
      # Estimating parameters for simulation
      mean_necropsies <- mean(mamData$Necropsies, na.rm = TRUE)
      variance_necropsies <- var(mamData$Necropsies, na.rm = TRUE)
      size_necropsy <- (mean_necropsies^2) / (variance_necropsies - mean_necropsies)
      prob_necropsy <- mean_necropsies / variance_necropsies
      
      # Simulating new necropsies data
      simulated_necropsies <- rnbinom(n, size = size_necropsy, mu = mean_necropsies)
      simulated_necropsies <- ifelse(simulated_necropsies < min_necropsies, min_necropsies, simulated_necropsies)
      mamData$SimulatedNecropsies <- simulated_necropsies
      



#adult weight models
#adult weight neo

mamData$SE_simpleSim <-1/sqrt(mamData$SimulatedNecropsies)

cutData <- mamData[,c(11,9,43,42,38),drop=FALSE] 
cutData$adult_weight[cutData$adult_weight < 0] <- NA
cutData <- na.omit(cutData)


cutData$Species <- gsub(" ", "_", cutData$Species) 
includedSpecies<-cutData$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species
SE<-setNames(cutData$SE_simpleSim,cutData$Species)[rownames(cutData)]

# Adjusting the simulation of 'simNeo' with increased variability
adult_weightFrac<-cutData$adult_weight/max(cutData$adult_weight)
noise_level <- sd(adult_weightFrac)*0.0000001 # Keep noise level consistent
intercept<-min(adult_weightFrac)
simNeo <- intercept+slope * adult_weightFrac + rnorm(nrow(cutData), mean = mean(adult_weightFrac), sd = noise_level)
simNeo[simNeo < 0] <- 0
cutData$simNeo<-simNeo
#pgls model
adult.weight.neo<-pglsSEyPagel(simNeo~log10(adult_weight),data=cutData,tree=pruned.tree,se=SE,method = "ML")

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



NeoplasiaOccurences<-round(as.numeric(cutData$simNeo*cutData$SimulatedNecropsies))
NonOccurences<-round(cutData$SimulatedNecropsies- NeoplasiaOccurences)
adultWeight<-log10(cutData$adult_weight)
pglsPagel_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "PGLSSEY", pval = p.v.adult.weight.neo,simslope= slope, slope = slope.adult.weight.neo)

tryCatch({
  compar <- suppressMessages(compar.gee(cbind(NeoplasiaOccurences, NonOccurences) ~ adultWeight, phy = pruned.tree, family = "binomial"))
  cNpval <- foo(compar)[2,4] # Assuming summary(compar) returns the correct structure
  cNslope <- foo(compar)[2,1]
  compar_row <- data.frame(n = nrow(cutData),min = min_necropsies, model = "Compar", pval = cNpval,simslope= slope,  slope = cNslope)
}, error = function(e) {
  # If an error occurs, assign NA to these variables
  compar_row <- data.frame(n = NA,min = NA, model = "Compar", pval = NA,simslope= NA,  slope = NA)
})
# ggplot(cutData, aes(x = adultWeight, y = NeoplasiaOccurences / (NeoplasiaOccurences + NonOccurences))) +
#   geom_point() +
#   labs(title = "Effect of Adult Weight on Neoplasia Occurrences", x = "Adult Weight", y = "Probability of Neoplasia Occurrence") +
#   theme_minimal()


all_results[[length(all_results) + 1]] <- rbind(pglsPagel_row, compar_row)

# Combine results for this simulation
#dataReport <- rbind(dataReport, pglsPagel_row, compar_row)

      }
    }
  }
  result_rows <- do.call(rbind, all_results)
  return(result_rows)
}




results_list <- mclapply(slope_values, run_simulation, mc.cores = no_cores)

# Combine all results from each slope into one data frame
final_results <- do.call(rbind, results_list)

pglsresults<-filter(dataReport, model == "PGLSSEY")
pglssig<-nrow(filter(pglsresults, pval <= 0.05))
levelpgls<-pglssig/nrow(pglsresults) 


gessesults<-filter(dataReport, model == "Compar")
gesssig<-nrow(filter(gessesults, pval <= 0.05))
levelgee<-gesssig/nrow(gessesults) 

mean(pglsresults$slope)
mean(gessesults$slope)


#sig levels for min 1 50 species
min1species<-filter(dataReport, min == 1 )
min1species50<-filter(min1species, n <=50 )
min150pglsr<-filter(min1species50, model == "PGLSSEY")
min1species50sigs<-nrow(filter(min150pglsr, pval <= 0.05))
min150levelpgls<-min1species50sigs/nrow(min150pglsr) 
min150geer<-filter(min1species50, model == "Compar")
min1species50sigsgee<-nrow(filter(min150geer, pval <= 0.05))
min150levelgee<-min1species50sigsgee/nrow(min150geer) 



#sig levels for min 1 75 species
min1species75<-filter(min1species, n >75, n <= 75)
min175pglsr<-filter(min1species75, model == "PGLSSEY")
min1species75sigs<-nrow(filter(min175pglsr, pval <= 0.05))
min175levelpgls<-min1species75sigs/nrow(min175pglsr) 
min175geer<-filter(min1species75, model == "Compar")
min1species75sigsgee<-nrow(filter(min175geer, pval <= 0.05))
min175levelgee<-min1species75sigsgee/nrow(min175geer) 

#sig levels for min 1 100 species
min1species100<-filter(min1species, n > 75, n <= 100)
min1100pglsr<-filter(min1species100, model == "PGLSSEY")
min1species100sigs<-nrow(filter(min1100pglsr, pval <= 0.05))
min1100levelpgls<-min1species100sigs/nrow(min1100pglsr) 
min1100pglsslope<-mean(min1100pglsr$slope)

min1100geer<-filter(min1species100, model == "Compar")
min1species100sigsgee<-nrow(filter(min1100geer, pval <= 0.05))
min1100levelgee<-min1species100sigsgee/nrow(min1100geer) 
min1100geeslope<-mean(min1100geer$slope)



#sig levels for min 1 200 species
min1species200<-filter(min1species, n > 100, n <= 200)
min1200pglsr<-filter(min1species200, model == "PGLSSEY")
min1species200sigs<-nrow(filter(min1200pglsr, pval <= 0.05))
min1200levelpgls<-min1species200sigs/nrow(min1200pglsr) 
min1200geer<-filter(min1species200, model == "Compar")
min1species200sigsgee<-nrow(filter(min1200geer, pval <= 0.05))
min1200levelgee<-min1species200sigsgee/nrow(min1200geer) 


#sig levels for min 1 300 species
min1species300<-filter(min1species, n > 200, n <= 300)
min1300pglsr<-filter(min1species300, model == "PGLSSEY")
min1species300sigs<-nrow(filter(min1300pglsr, pval <= 0.05))
min1300levelpgls<-min1species300sigs/nrow(min1300pglsr) 
min1300geer<-filter(min1species300, model == "Compar")
min1species300sigsgee<-nrow(filter(min1300geer, pval <= 0.05))
min1300levelgee<-min1species300sigsgee/nrow(min1300geer) 

#sig levels for min 1 400 species
min1species400<-filter(min1species, n > 300, n <= 400)
min1400pglsr<-filter(min1species400, model == "PGLSSEY")
min1species400sigs<-nrow(filter(min1400pglsr, pval <= 0.05))
min1400levelpgls<-min1species400sigs/nrow(min1400pglsr) 
min1400geer<-filter(min1species400, model == "Compar")
min1species400sigsgee<-nrow(filter(min1400geer, pval <= 0.05))
min1400levelgee<-min1species400sigsgee/nrow(min1400geer) 





#sig levels for min 2 50 species
min20species<-filter(dataReport, min == 20 )
min20species50<-filter(min20species, n <=50 )
min2050pglsr<-filter(min20species50, model == "PGLSSEY")
min20species50sigs<-nrow(filter(min2050pglsr, pval <= 0.05))
min2050levelpgls<-min20species50sigs/nrow(min2050pglsr) 
min2050geer<-filter(min20species50, model == "Compar")
min20species50sigsgee<-nrow(filter(min2050geer, pval <= 0.05))
min2050levelgee<-min20species50sigsgee/nrow(min2050geer) 


#sig levels for min 1 75 species
min20species75<-filter(min20species, n >50, n <= 75)
min2075pglsr<-filter(min20species75, model == "PGLSSEY")
min20species75sigs<-nrow(filter(min2075pglsr, pval <= 0.05))
min2075levelpgls<-min20species75sigs/nrow(min2075pglsr) 
min2075geer<-filter(min20species75, model == "Compar")
min20species75sigsgee<-nrow(filter(min2075geer, pval <= 0.05))
min2075levelgee<-min20species75sigsgee/nrow(min2075geer) 

#sig levels for min 20 100 species
min20species100<-filter(min20species, n > 75, n <= 100)
min20100pglsr<-filter(min20species100, model == "PGLSSEY")
min20species100sigs<-nrow(filter(min20100pglsr, pval <= 0.05))
min20100levelpgls<-min20species100sigs/nrow(min20100pglsr) 
min20100geer<-filter(min20species100, model == "Compar")
min20species100sigsgee<-nrow(filter(min20100geer, pval <= 0.05))
min20100levelgee<-min20species100sigsgee/nrow(min20100geer) 


#sig levels for min 1 200 species
min20species200<-filter(min20species, n > 100, n <= 200)
min20200pglsr<-filter(min20species200, model == "PGLSSEY")
min20species200sigs<-nrow(filter(min20200pglsr, pval <= 0.05))
min20200levelpgls<-min20species200sigs/nrow(min20200pglsr) 
min20200geer<-filter(min20species200, model == "Compar")
min20species200sigsgee<-nrow(filter(min20200geer, pval <= 0.05))
min20200levelgee<-min20species200sigsgee/nrow(min20200geer) 


#sig levels for min 1 300 species
min20species300<-filter(min20species, n > 200, n <= 300)
min20300pglsr<-filter(min20species300, model == "PGLSSEY")
min20species300sigs<-nrow(filter(min20300pglsr, pval <= 0.05))
min20300levelpgls<-min20species300sigs/nrow(min20300pglsr) 
min20300geer<-filter(min20species300, model == "Compar")
min20species300sigsgee<-nrow(filter(min20300geer, pval <= 0.05))
min20300levelgee<-min20species300sigsgee/nrow(min20300geer) 

#sig levels for min 1 400 species
min20species400<-filter(min20species, n > 300, n <= 400)
min20400pglsr<-filter(min20species400, model == "PGLSSEY")
min20species400sigs<-nrow(filter(min20400pglsr, pval <= 0.05))
min20400levelpgls<-min20species400sigs/nrow(min20400pglsr) 
min20400geer<-filter(min20species400, model == "Compar")
min20species400sigsgee<-nrow(filter(min20400geer, pval <= 0.05))
min20400levelgee<-min20species400sigsgee/nrow(min20400geer) 



write.csv(dataReport, file = "./expExperimentsSols.csv")

