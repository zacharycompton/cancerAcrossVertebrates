#open libraries
library(tidyverse)
library(phytools)

#read tree
tree <- read.tree("min20Fixed516.nwk")
tree
plot(tree)

## plot tree
plotTree(tree,type="fan", ftype = "i")

## read csv file for mammals

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Mammalia")))
row.names(Data)= gsub(" ", "_", row.names(Data))

plot(tree)

## prune the tree to match the data
includedSpecies <- row.names(Data)
newtips<-str_remove_all(tree$tip.label,"_ott")
newtips<-str_remove_all(newtips,".ott")
newtips<-str_remove_all(newtips,"-ott")
newtips<-str_remove_all(newtips,"[1234567890]")
newtips<-sub('^([^_]+_[^_]+).*', '\\1', newtips)

## pruning the tree
tree$tip.label <- newtips
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)


plot(pruned.tree)

## removing discrepancies
Data$Keep <- row.names(Data) %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


plotTree(pruned.tree,type="fan",ftype="i")

## match species from csv to the tree
neoplasia=dplyr::filter(Data, row.names(Data) %in% pruned.tree$tip.label)

## select trait
neoplasia=neoplasia%>%dplyr::select(NeoplasiaPrevalence)

## change this into a vector
neoplasia=as.matrix(neoplasia)[,1]
neoplasia

## estimate ancestral states
fit<-fastAnc(pruned.tree,neoplasia,vars=TRUE,CI=TRUE)
fit

## confidence intervals
fit$CI[1,]

## projection of the reconstruction onto the edges of the tree
Mobj<-contMap(pruned.tree,neoplasia,plot=TRUE,res=200, type = "fan")

## invert color scale so red is highest trait value and blue is lowest
Mobj<-setMap(Mobj,invert=TRUE)

## plot 
m<-plot(Mobj,
     fsize=c(0.4,0.8),leg.txt="neoplasia prevalence")





## read csv file for Sauropsids

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
row.names(Data)= gsub(" ", "_", row.names(Data))

plot(tree)

## prune the tree to match the data
includedSpecies <- row.names(Data)
newtips<-str_remove_all(tree$tip.label,"_ott")
newtips<-str_remove_all(newtips,".ott")
newtips<-str_remove_all(newtips,"-ott")
newtips<-str_remove_all(newtips,"[1234567890]")
newtips<-sub('^([^_]+_[^_]+).*', '\\1', newtips)

## pruning the tree
tree$tip.label <- newtips
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)


plot(pruned.tree)

## removing discrepancies
Data$Keep <- row.names(Data) %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


plotTree(pruned.tree,type="fan",ftype="i")

## match species from csv to the tree
neoplasia=dplyr::filter(Data, row.names(Data) %in% pruned.tree$tip.label)

## select trait
neoplasia=neoplasia%>%dplyr::select(NeoplasiaPrevalence)

## change this into a vector
neoplasia=as.matrix(neoplasia)[,1]
neoplasia

## estimate ancestral states
fit<-fastAnc(pruned.tree,neoplasia,vars=TRUE,CI=TRUE)
fit

## confidence intervals
fit$CI[1,]

## projection of the reconstruction onto the edges of the tree
Sobj<-contMap(pruned.tree,neoplasia,plot=TRUE,res=200, type = "fan")

## invert color scale so red is highest trait value and blue is lowest
Sobj<-setMap(Sobj,invert=TRUE)

## plot 
s<-plot(Sobj,
           fsize=c(0.4,0.8),leg.txt="neoplasia prevalence")





## read csv file for Amphibians

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Amphibia")))
row.names(Data)= gsub(" ", "_", row.names(Data))

plot(tree)

## prune the tree to match the data
includedSpecies <- row.names(Data)
newtips<-str_remove_all(tree$tip.label,"_ott")
newtips<-str_remove_all(newtips,".ott")
newtips<-str_remove_all(newtips,"-ott")
newtips<-str_remove_all(newtips,"[1234567890]")
newtips<-sub('^([^_]+_[^_]+).*', '\\1', newtips)

## pruning the tree
tree$tip.label <- newtips
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)


plot(pruned.tree)

## removing discrepancies
Data$Keep <- row.names(Data) %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


plotTree(pruned.tree,type="fan",ftype="i")

## match species from csv to the tree
neoplasia=dplyr::filter(Data, row.names(Data) %in% pruned.tree$tip.label)

## select trait
neoplasia=neoplasia%>%dplyr::select(NeoplasiaPrevalence)

## change this into a vector
neoplasia=as.matrix(neoplasia)[,1]
neoplasia

## estimate ancestral states
fit<-fastAnc(pruned.tree,neoplasia,vars=TRUE,CI=TRUE)
fit

## confidence intervals
fit$CI[1,]

## projection of the reconstruction onto the edges of the tree
Aobj<-contMap(pruned.tree,neoplasia,plot=TRUE,res=200, type = "fan")

## invert color scale so red is highest trait value and blue is lowest
Aobj<-setMap(Aobj,invert=TRUE)

## plot multiple phylos on to one image

layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
par(mar=c(4, 4, 4, 4))
plot(Sobj,
     fsize=c(0.3,0.8),leg.txt="neoplasia prevalence")
title(main="Sauropsids", cex.main=1, line = -.5)
par(mar=c(4, 4, 4, 4))
plot(Mobj,
     fsize=c(0.4,0.8),leg.txt="neoplasia prevalence")
title(main="Mammals", cex.main=1, line = -.5)
par(mar=c(4, 4, 4, 4))
plot(Aobj,
     fsize=c(0.4,0.8),leg.txt="neoplasia prevalence")
title(main="Amphibians", cex.main=1, line = -.5)




