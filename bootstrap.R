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
#library(xlsx)




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


Data <- read.csv("min20-2022.05.16.csv")
cleanPath<-read.csv("allrecords.csv")
tree <- read.tree("min20Fixed516.nwk")

cleanPath <- cleanPath %>% filter(Necropsy != -1)
cleanPath <- cleanPath %>% filter(Necropsy != 0)

# Convert both lists to lowercase for case-insensitive comparison
lowerInd <- tolower(cleanPath$Species)
lowerSpecies <- tolower(Data$Species)

# Keep only records that match species names
cleanPath <- cleanPath[lowerInd %in% lowerSpecies,]


distinct_species_count <- cleanPath %>%
  select(Species) %>%
  distinct() %>%
  nrow()


species_counts <- cleanPath %>%
  count(Species)


# Set the number of repetitions
num_repetitions <- 100

# Loop through each species
unique_species <- unique(cleanPath$Species)


valuesWeightMal<-data.frame(pvalues = NULL, R = NULL, Lambda = NULL)
valuesWeightNeo<-data.frame(pvalues = NULL, R = NULL, Lambda = NULL)

for (i in 1:num_repetitions) {
  # Create an empty data frame to store the results
  sampled_data <- data.frame()
  for (species_name in unique_species) {
    # Subset the data for the current species
    species_subset <- cleanPath %>%
      filter(Species == species_name)
    
    # For each individual ID, select the record with the highest value for Malignant
    highest_malignant_per_individual <- species_subset %>%
      group_by(ID) %>%
      slice(which.max(Malignant)) %>%
      ungroup()
    
    # Randomly sample half of the records for the current species
    half_sample <- highest_malignant_per_individual %>%
      sample_n(nrow(highest_malignant_per_individual) / 2)
    
    # Append the sampled data to the results data frame
    sampled_data <- bind_rows(sampled_data, half_sample)
  }
  # Reset row names in the resulting data frame
  rownames(sampled_data) <- NULL
  
  currentData <- sampled_data[,c(18,13,19,28,46),drop=FALSE]
  currentData$Type[is.na(currentData$Type)] <- "none"
  
  
  species_summary <- currentData %>%
    group_by(Species) %>%
    summarise(
      total_records = n(),
      sum_Malignant = sum(Malignant[Type != "cyst"]),
      sum_Masspresent = n()-sum(Type == "none")-sum(Type == "cyst"), 
      max_adult_weight = max(adult_weight) # Replace 'adult_weight' with the actual column name
    ) %>%
    ungroup()
  
  species_summary <- species_summary %>%
    mutate(SE_Simple = 1 / sqrt(total_records))
  
  species_summary <- species_summary %>%
    mutate(MalignancyPrevalence = sum_Malignant/total_records)
  
  species_summary <- species_summary %>%
    mutate(NeoplasiaPrevalence = sum_Masspresent/total_records)
  
  # Create a data frame with old and new species names

    old_species = c(
      "Ailurus_fulgens", "Giraffa_camelopardalis", "Chinchilla_chinchilla",
      "Egretta_garzetta", "Ciconia_abdimii", "Phoenicopterus_chilensis",
      "Phoenicopterus_ruber", "Athene_cunicularia"
    )
    new_species = c(
      "Ailurus", "Giraffa_reticulata", "Chinchilla_lanigera",
      "Egretta_novaehollandiae", "Ciconia_ciconia", "Phoenicopterus",
      "Phoenicopterus", "Athene")
  
  for (nameSwitch in 1: length(new_species)){
    species_summary <- species_summary %>%
      mutate(Species = if_else(Species == old_species[nameSwitch], new_species[nameSwitch], Species))
    
    
  }
  
    species_summary <- species_summary %>%
      mutate(max_adult_weight = ifelse(max_adult_weight == -1, NA, max_adult_weight))
    species_summary<-na.omit(species_summary)
  
  species_summary$Species <- gsub(" ", "_", species_summary$Species) 
  includedSpecies<-species_summary$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  species_summary$Keep <- species_summary$Species %in% pruned.tree$tip.label
  species_summary <- species_summary[!(species_summary$Keep==FALSE),]
  rownames(species_summary)<-species_summary$Species
  SE<-setNames(species_summary$SE_Simple,species_summary$Species)[rownames(species_summary)]
  
  weightPGLSmal<-pglsSEyPagel(MalignancyPrevalence~log10(max_adult_weight),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simwgt.mal <- summary(weightPGLSmal)$corBeta
  r.v.simwgt.mal <- format(r.v.simwgt.mal[2,1])
  r.v.simwgt.mal <-signif(as.numeric(r.v.simwgt.mal)^2, digits= 2)
  ld.v.simwgt.mal<- summary(weightPGLSmal)$modelStruct$corStruct
  ld.v.simwgt.mal <- signif(ld.v.simwgt.mal[1], digits = 2)
  p.v.simwgt.mal<-summary(weightPGLSmal)$tTable
  p.v.simwgt.mal<-signif(p.v.simwgt.mal[2,4], digits = 2)
  
  
  weightPGLSneo<-pglsSEyPagel(NeoplasiaPrevalence~log10(max_adult_weight),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simwgt.neo <- summary(weightPGLSneo)$corBeta
  r.v.simwgt.neo <- format(r.v.simwgt.neo[2,1])
  r.v.simwgt.neo <-signif(as.numeric(r.v.simwgt.neo)^2, digits= 2)
  ld.v.simwgt.neo<- summary(weightPGLSneo)$modelStruct$corStruct
  ld.v.simwgt.neo <- signif(ld.v.simwgt.neo[1], digits = 2)
  p.v.simwgt.neo<-summary(weightPGLSneo)$tTable
  p.v.simwgt.neo<-signif(p.v.simwgt.neo[2,4], digits = 2)
  
  
  # Create a new row with a value to add at the end
  new_row <- data.frame(pvalues = p.v.simwgt.mal, R = r.v.simwgt.mal, Lambda = ld.v.simwgt.mal)
  
  # Add the new row at the end of the dataframe
  valuesWeightMal <- rbind(valuesWeightMal, new_row)
  
  
  # Create a new row with a value to add at the end
  new_rowNeo <- data.frame(pvalues = p.v.simwgt.neo, R = r.v.simwgt.neo, Lambda = ld.v.simwgt.neo)
  
  # Add the new row at the end of the dataframe
  valuesWeightNeo <- rbind(valuesWeightNeo, new_rowNeo)
  
  
}


#Gestation



Data <- read.csv("min20-2022.05.16.csv")
cleanPath<-read.csv("allrecords.csv")
tree <- read.tree("min20Fixed516.nwk")

cleanPath <- cleanPath %>% filter(Necropsy != -1)
cleanPath <- cleanPath %>% filter(Necropsy != 0)

 #Convert both lists to lowercase for case-insensitive comparison
lowerInd <- tolower(cleanPath$Species)
lowerSpecies <- tolower(Data$Species)

# Keep only records that match species names
cleanPath <- cleanPath[lowerInd %in% lowerSpecies,]


distinct_species_count <- cleanPath %>%
  select(Species) %>%
  distinct() %>%
  nrow()


species_counts <- cleanPath %>%
  count(Species)


# Set the number of repetitions
num_repetitions <- 100

# Loop through each species
unique_species <- unique(cleanPath$Species)


valuesGestMal<-data.frame(pvalues = NULL, R = NULL, Lambda = NULL)
valuesGestNeo<-data.frame(pvalues = NULL, R = NULL, Lambda = NULL)

for (i in 1:num_repetitions) {
  # Create an empty data frame to store the results
  sampled_data <- data.frame()
  for (species_name in unique_species) {
    # Subset the data for the current species
    species_subset <- cleanPath %>%
      filter(Species == species_name)
    
    # For each individual ID, select the record with the highest value for Malignant
    highest_malignant_per_individual <- species_subset %>%
      group_by(ID) %>%
      slice(which.max(Malignant)) %>%
      ungroup()
    
    # Randomly sample half of the records for the current species
    half_sample <- highest_malignant_per_individual %>%
      sample_n(nrow(highest_malignant_per_individual) / 2)
    
    # Append the sampled data to the results data frame
    sampled_data <- bind_rows(sampled_data, half_sample)
  }
  # Reset row names in the resulting data frame
  rownames(sampled_data) <- NULL
  
  sampled_data$Malignant[sampled_data$Malignant == -1] <- 0
  
  currentData <- sampled_data[,c(18,13,19,28,38),drop=FALSE]
  currentData$Type[is.na(currentData$Type)] <- "none"
  
  
  
  #currentData <- currentData %>%
    #mutate_all(~ifelse(. == -1, NA, .))
  
  #currentData<-na.omit(currentData)
  
  
  species_summary <- currentData %>%
    group_by(Species) %>%
    summarise(
      total_records = n(),
      sum_Malignant = sum(Malignant[Type != "cyst"]),
      sum_Masspresent = n()-sum(Type == "none")-sum(Type == "cyst"),  
      max_Gestation = max(Gestation)
    ) %>%
    ungroup()
  
  species_summary <- species_summary %>%
    mutate(SE_Simple = 1 / sqrt(total_records))
  
  species_summary <- species_summary %>%
    mutate(MalignancyPrevalence = sum_Malignant/total_records)
  
  species_summary <- species_summary %>%
    mutate(NeoplasiaPrevalence = sum_Masspresent/total_records)
  
  
  # Create a data frame with old and new species names
  
  old_species = c(
    "Ailurus_fulgens", "Giraffa_camelopardalis", "Chinchilla_chinchilla",
    "Egretta_garzetta", "Ciconia_abdimii", "Phoenicopterus_chilensis",
    "Phoenicopterus_ruber", "Athene_cunicularia"
  )
  new_species = c(
    "Ailurus", "Giraffa_reticulata", "Chinchilla_lanigera",
    "Egretta_novaehollandiae", "Ciconia_ciconia", "Phoenicopterus",
    "Phoenicopterus", "Athene")
  
  for (nameSwitch in 1: length(new_species)){
    species_summary <- species_summary %>%
      mutate(Species = if_else(Species == old_species[nameSwitch], new_species[nameSwitch], Species))
    
    
  }
  
  
  
  
  species_summary <- species_summary %>%
    mutate(max_Gestation = ifelse(max_Gestation == -1, NA, max_Gestation))
  species_summary<-na.omit(species_summary)
  
  species_summary$Species <- gsub(" ", "_", species_summary$Species) 
  includedSpecies<-species_summary$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  species_summary$Keep <- species_summary$Species %in% pruned.tree$tip.label
  species_summary <- species_summary[!(species_summary$Keep==FALSE),]
  rownames(species_summary)<-species_summary$Species
  SE<-setNames(species_summary$SE_Simple,species_summary$Species)[rownames(species_summary)]
  species_summary$max_Gestation
  
  
  GestPGLSmal<-pglsSEyPagel(MalignancyPrevalence~log10(max_Gestation),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simgest.mal <- summary(GestPGLSmal)$corBeta
  r.v.simgest.mal <- format(r.v.simgest.mal[2,1])
  r.v.simgest.mal <-signif(as.numeric(r.v.simgest.mal)^2, digits= 2)
  ld.v.simgest.mal<- summary(GestPGLSmal)$modelStruct$corStruct
  ld.v.simgest.mal <- signif(ld.v.simgest.mal[1], digits = 2)
  p.v.simgest.mal<-summary(GestPGLSmal)$tTable
  p.v.simgest.mal<-signif(p.v.simgest.mal[2,4], digits = 2)
  
  
  GestPGLSneo<-pglsSEyPagel(NeoplasiaPrevalence~log10(max_Gestation),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simgest.neo <- summary(GestPGLSneo)$corBeta
  r.v.simgest.neo <- format(r.v.simgest.neo[2,1])
  r.v.simgest.neo <-signif(as.numeric(r.v.simgest.neo)^2, digits= 2)
  ld.v.simgest.neo<- summary(GestPGLSneo)$modelStruct$corStruct
  ld.v.simgest.neo <- signif(ld.v.simgest.neo[1], digits = 2)
  p.v.simgest.neo<-summary(GestPGLSneo)$tTable
  p.v.simgest.neo<-signif(p.v.simgest.neo[2,4], digits = 2)
  
  
  # Create a new row with a value to add at the end
  new_row <- data.frame(pvalues = p.v.simgest.mal, R = r.v.simgest.mal, Lambda = ld.v.simgest.mal)
  
  # Add the new row at the end of the dataframe
  valuesGestMal <- rbind(valuesGestMal, new_row)
  
  
  # Create a new row with a value to add at the end
  new_rowNeo <- data.frame(pvalues = p.v.simgest.neo, R = r.v.simgest.neo, Lambda = ld.v.simgest.neo)
  
  # Add the new row at the end of the dataframe
  valuesGestNeo <- rbind(valuesGestNeo, new_rowNeo)
  
  
}




#longevity



Data <- read.csv("min20-2022.05.16.csv")
cleanPath<-read.csv("allrecords.csv")
tree <- read.tree("min20Fixed516.nwk")


cleanPath <- cleanPath %>% filter(Necropsy != -1)
cleanPath <- cleanPath %>% filter(Necropsy != 0)


# Convert both lists to lowercase for case-insensitive comparison
lowerInd <- tolower(cleanPath$Species)
lowerSpecies <- tolower(Data$Species)

# Keep only records that match species names
cleanPath <- cleanPath[lowerInd %in% lowerSpecies,]


distinct_species_count <- cleanPath %>%
  select(Species) %>%
  distinct() %>%
  nrow()


species_counts <- cleanPath %>%
  count(Species)


# Set the number of repetitions
num_repetitions <- 100

# Loop through each species
unique_species <- unique(cleanPath$Species)


valuesLongMal<-data.frame(pvalues = NULL, R = NULL, Lambda = NULL)
valuesLongNeo<-data.frame(pvalues = NULL, R = NULL, Lambda = NULL)

for (i in 1:num_repetitions) {
  # Create an empty data frame to store the results
  sampled_data <- data.frame()
  for (species_name in unique_species) {
    # Subset the data for the current species
    species_subset <- cleanPath %>%
      filter(Species == species_name)
    
    # For each individual ID, select the record with the highest value for Malignant
    highest_malignant_per_individual <- species_subset %>%
      group_by(ID) %>%
      slice(which.max(Malignant)) %>%
      ungroup()
    
    # Randomly sample half of the records for the current species
    half_sample <- highest_malignant_per_individual %>%
      sample_n(nrow(highest_malignant_per_individual) / 2)
    
    # Append the sampled data to the results data frame
    sampled_data <- bind_rows(sampled_data, half_sample)
  }
  # Reset row names in the resulting data frame
  rownames(sampled_data) <- NULL
  
  currentData <- sampled_data[,c(18,13,19,28,48),drop=FALSE]
  
  
  currentData <- currentData %>%
    mutate_all(~ifelse(. == -1, NA, .))
  
  currentData<-na.omit(currentData)
  
  
  species_summary <- currentData %>%
    group_by(Species) %>%
    summarise(
      total_records = n(),
      sum_Malignant = sum(Malignant[Type != "cyst"]),
      sum_Masspresent = n()-sum(Type == "none")-sum(Type == "cyst"),  
      max_longevity = max(max_longevity) # Replace 'adult_weight' with the actual column name
    ) %>%
    ungroup()
  
  species_summary <- species_summary %>%
    mutate(SE_Simple = 1 / sqrt(total_records))
  
  species_summary <- species_summary %>%
    mutate(MalignancyPrevalence = sum_Malignant/total_records)
  
  species_summary <- species_summary %>%
    mutate(NeoplasiaPrevalence = sum_Masspresent/total_records)
  
  # Create a data frame with old and new species names
  
  old_species = c(
    "Ailurus_fulgens", "Giraffa_camelopardalis", "Chinchilla_chinchilla",
    "Egretta_garzetta", "Ciconia_abdimii", "Phoenicopterus_chilensis",
    "Phoenicopterus_ruber", "Athene_cunicularia"
  )
  new_species = c(
    "Ailurus", "Giraffa_reticulata", "Chinchilla_lanigera",
    "Egretta_novaehollandiae", "Ciconia_ciconia", "Phoenicopterus",
    "Phoenicopterus", "Athene")
  
  for (nameSwitch in 1: length(new_species)){
    species_summary <- species_summary %>%
      mutate(Species = if_else(Species == old_species[nameSwitch], new_species[nameSwitch], Species))
    
    
  }
  
  species_summary <- species_summary %>%
    mutate(max_longevity = ifelse(max_longevity == -1, NA, max_longevity))
  species_summary<-na.omit(species_summary)
  
  species_summary$Species <- gsub(" ", "_", species_summary$Species) 
  includedSpecies<-species_summary$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  species_summary$Keep <- species_summary$Species %in% pruned.tree$tip.label
  species_summary <- species_summary[!(species_summary$Keep==FALSE),]
  rownames(species_summary)<-species_summary$Species
  SE<-setNames(species_summary$SE_Simple,species_summary$Species)[rownames(species_summary)]
  
  LongPGLSmal<-pglsSEyPagel(MalignancyPrevalence~log10(max_longevity),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simlong.mal <- summary(LongPGLSmal)$corBeta
  r.v.simlong.mal <- format(r.v.simlong.mal[2,1])
  r.v.simlong.mal <-signif(as.numeric(r.v.simlong.mal)^2, digits= 2)
  ld.v.simlong.mal<- summary(LongPGLSmal)$modelStruct$corStruct
  ld.v.simlong.mal <- signif(ld.v.simlong.mal[1], digits = 2)
  p.v.simlong.mal<-summary(LongPGLSmal)$tTable
  p.v.simlong.mal<-signif(p.v.simlong.mal[2,4], digits = 2)
  
  
  LongPGLSneo<-pglsSEyPagel(NeoplasiaPrevalence~log10(max_longevity),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simlong.neo <- summary(LongPGLSneo)$corBeta
  r.v.simlong.neo <- format(r.v.simlong.neo[2,1])
  r.v.simlong.neo <-signif(as.numeric(r.v.simlong.neo)^2, digits= 2)
  ld.v.simlong.neo<- summary(LongPGLSneo)$modelStruct$corStruct
  ld.v.simlong.neo <- signif(ld.v.simlong.neo[1], digits = 2)
  p.v.simlong.neo<-summary(LongPGLSneo)$tTable
  p.v.simlong.neo<-signif(p.v.simlong.neo[2,4], digits = 2)
  
  
  # Create a new row with a value to add at the end
  new_row <- data.frame(pvalues = p.v.simlong.mal, R = r.v.simlong.mal, Lambda = ld.v.simlong.mal)
  
  # Add the new row at the end of the dataframe
  valuesLongMal <- rbind(valuesLongMal, new_row)
  
  
  # Create a new row with a value to add at the end
  new_rowNeo <- data.frame(pvalues = p.v.simlong.neo, R = r.v.simlong.neo, Lambda = ld.v.simlong.neo)
  
  # Add the new row at the end of the dataframe
  valuesLongNeo <- rbind(valuesLongNeo, new_rowNeo)
  
  
}



#Gestation + BodyMass



Data <- read.csv("min20-2022.05.16.csv")
cleanPath<-read.csv("allrecords.csv")
tree <- read.tree("min20Fixed516.nwk")

cleanPath <- cleanPath %>% filter(Necropsy != -1)
cleanPath <- cleanPath %>% filter(Necropsy != 0)

#Convert both lists to lowercase for case-insensitive comparison
lowerInd <- tolower(cleanPath$Species)
lowerSpecies <- tolower(Data$Species)

# Keep only records that match species names
cleanPath <- cleanPath[lowerInd %in% lowerSpecies,]


distinct_species_count <- cleanPath %>%
  select(Species) %>%
  distinct() %>%
  nrow()


species_counts <- cleanPath %>%
  count(Species)


# Set the number of repetitions
num_repetitions <- 100

# Loop through each species
unique_species <- unique(cleanPath$Species)


valuesGestMassMal<-data.frame(pvalGest = NULL,pvalMass = NULL, R = NULL, Lambda = NULL)
valuesGestMassNeo<-data.frame(pvalGest = NULL,pvalMass = NULL, R = NULL, Lambda = NULL)

for (i in 1:num_repetitions) {
  # Create an empty data frame to store the results
  sampled_data <- data.frame()
  for (species_name in unique_species) {
    # Subset the data for the current species
    species_subset <- cleanPath %>%
      filter(Species == species_name)
    
    # For each individual ID, select the record with the highest value for Malignant
    highest_malignant_per_individual <- species_subset %>%
      group_by(ID) %>%
      slice(which.max(Malignant)) %>%
      ungroup()
    
    # Randomly sample half of the records for the current species
    half_sample <- highest_malignant_per_individual %>%
      sample_n(nrow(highest_malignant_per_individual) / 2)
    
    # Append the sampled data to the results data frame
    sampled_data <- bind_rows(sampled_data, half_sample)
  }
  # Reset row names in the resulting data frame
  rownames(sampled_data) <- NULL
  
  sampled_data$Malignant[sampled_data$Malignant == -1] <- 0
  
  currentData <- sampled_data[,c(18,13,19,28,46,38),drop=FALSE]
  currentData$Type[is.na(currentData$Type)] <- "none"
  
  
  
  #currentData <- currentData %>%
  #mutate_all(~ifelse(. == -1, NA, .))
  
  #currentData<-na.omit(currentData)
  
  
  species_summary <- currentData %>%
    group_by(Species) %>%
    summarise(
      total_records = n(),
      sum_Malignant = sum(Malignant[Type != "cyst"]),
      sum_Masspresent = n()-sum(Type == "none")-sum(Type == "cyst"),  
      max_Gestation = max(Gestation),
      max_Mass = max(adult_weight)
    ) %>%
    ungroup()
  
  species_summary <- species_summary %>%
    mutate(SE_Simple = 1 / sqrt(total_records))
  
  species_summary <- species_summary %>%
    mutate(MalignancyPrevalence = sum_Malignant/total_records)
  
  species_summary <- species_summary %>%
    mutate(NeoplasiaPrevalence = sum_Masspresent/total_records)
  
  
  # Create a data frame with old and new species names
  
  old_species = c(
    "Ailurus_fulgens", "Giraffa_camelopardalis", "Chinchilla_chinchilla",
    "Egretta_garzetta", "Ciconia_abdimii", "Phoenicopterus_chilensis",
    "Phoenicopterus_ruber", "Athene_cunicularia"
  )
  new_species = c(
    "Ailurus", "Giraffa_reticulata", "Chinchilla_lanigera",
    "Egretta_novaehollandiae", "Ciconia_ciconia", "Phoenicopterus",
    "Phoenicopterus", "Athene")
  
  for (nameSwitch in 1: length(new_species)){
    species_summary <- species_summary %>%
      mutate(Species = if_else(Species == old_species[nameSwitch], new_species[nameSwitch], Species))
    
    
  }
  
  
  
  
  species_summary <- species_summary %>%
    mutate(max_Gestation = ifelse(max_Gestation == -1, NA, max_Gestation))%>%
    mutate(max_Mass = ifelse(max_Mass == -1, NA, max_Mass))
  species_summary<-na.omit(species_summary)
  
  species_summary$Species <- gsub(" ", "_", species_summary$Species) 
  includedSpecies<-species_summary$Species
  pruned.tree<-drop.tip(
    tree, setdiff(
      tree$tip.label, includedSpecies))
  pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
  species_summary$Keep <- species_summary$Species %in% pruned.tree$tip.label
  species_summary <- species_summary[!(species_summary$Keep==FALSE),]
  rownames(species_summary)<-species_summary$Species
  SE<-setNames(species_summary$SE_Simple,species_summary$Species)[rownames(species_summary)]
  
  
  GestMassPGLSmal<-pglsSEyPagel(MalignancyPrevalence~log10(max_Gestation)+log10(max_Mass),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simgestmass.mal <- summary(GestMassPGLSmal)$corBeta
  r.v.simgestmass.mal <- format(r.v.simgestmass.mal[2,1])
  r.v.simgestmass.mal <-signif(as.numeric(r.v.simgestmass.mal)^2, digits= 2)
  ld.v.simgestmass.mal<- summary(GestMassPGLSmal)$modelStruct$corStruct
  ld.v.simgestmass.mal <- signif(ld.v.simgestmass.mal[1], digits = 2)
  p.v.simgest.mal<-summary(GestMassPGLSmal)$tTable
  p.v.simgest.mal<-signif(p.v.simgest.mal[2,4], digits = 2)
  p.v.simmass.mal<-summary(GestMassPGLSmal)$tTable
  p.v.simmass.mal<-signif(p.v.simmass.mal[3,4], digits = 2)
  
  
  GestMassPGLSneo<-pglsSEyPagel(NeoplasiaPrevalence~log10(max_Gestation)+log10(max_Mass),data=species_summary,tree=pruned.tree,se=SE,method = "ML")
  
  r.v.simgestmass.neo <- summary(GestMassPGLSneo)$corBeta
  r.v.simgestmass.neo <- format(r.v.simgestmass.neo[2,1])
  r.v.simgestmass.neo <-signif(as.numeric(r.v.simgestmass.neo)^2, digits= 2)
  ld.v.simgestmass.neo<- summary(GestMassPGLSneo)$modelStruct$corStruct
  ld.v.simgestmass.neo <- signif(ld.v.simgestmass.neo[1], digits = 2)
  p.v.simgest.neo<-summary(GestMassPGLSneo)$tTable
  p.v.simgest.neo<-signif(p.v.simgest.neo[2,4], digits = 2)
  p.v.simmass.neo<-summary(GestMassPGLSneo)$tTable
  p.v.simmass.neo<-signif(p.v.simmass.neo[3,4], digits = 2)
  
  

  
  # Create a new row with a value to add at the end
  new_row <- data.frame(pvalGest = p.v.simgest.mal, pvalMass = p.v.simmass.mal,R = r.v.simgestmass.mal, Lambda = ld.v.simgestmass.mal)
  
  # Add the new row at the end of the dataframe
  valuesGestMassMal <- rbind(valuesGestMassMal, new_row)
  
  
  # Create a new row with a value to add at the end
  new_rowNeo <- data.frame(pvalGest = p.v.simgest.neo, pvalMass = p.v.simmass.neo,R = r.v.simgestmass.neo, Lambda = ld.v.simgestmass.neo)
  
  # Add the new row at the end of the dataframe
  valuesGestMassNeo <- rbind(valuesGestMassNeo, new_rowNeo)
  
  
}

