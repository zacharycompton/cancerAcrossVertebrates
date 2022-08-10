#open libraries
library(tidyverse)
library(phytools)
library(patchwork)
library(tidyverse)
library(cowplot)



par(mfrow=c(3,1))
#read tree
tree <- read.tree("min20Fixed516.nwk")

## read csv file for mammals

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Mammalia")))
row.names(Data)= gsub(" ", "_", row.names(Data))


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


## removing discrepancies
Data$Keep <- row.names(Data) %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


## match species from csv to the tree
neoplasia=dplyr::filter(Data, row.names(Data) %in% pruned.tree$tip.label)

## select trait
neoplasia=neoplasia%>%dplyr::select(NeoplasiaPrevalence)

## change this into a vector
neoplasia=as.matrix(neoplasia)[,1]
neoplasia

## estimate ancestral states
fit<-fastAnc(pruned.tree,neoplasia,vars=TRUE,CI=TRUE)

## confidence intervals
fit$CI[1,]

## projection of the reconstruction onto the edges of the tree
Mobj<-contMap(pruned.tree,neoplasia,plot=FALSE,res=300, type = "phylogram",lwd = 2)


## invert color scale so red is highest trait value and blue is lowest
Mobj<-setMap(Mobj,invert=TRUE)


## plot 
m<-plot(Mobj,
     fsize=c(1,0.8),leg.txt="neoplasia prevalence",lwd= 2.5)
title(main="A", adj = .05, line = -1)



## read csv file for Sauropsids

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
row.names(Data)= gsub(" ", "_", row.names(Data))


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


## removing discrepancies
Data$Keep <- row.names(Data) %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


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
Sobj<-contMap(pruned.tree,neoplasia,plot=FALSE,res=200, type = "fan")

## invert color scale so red is highest trait value and blue is lowest
Sobj<-setMap(Sobj,invert=TRUE)

## plot 
s<-plot(Sobj,
           fsize=c(1,0.8),leg.txt="neoplasia prevalence",lwd = 2.5)
title(main="B", adj = .05, line = -1)





## read csv file for Amphibians

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Amphibia")))
row.names(Data)= gsub(" ", "_", row.names(Data))


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


## removing discrepancies
Data$Keep <- row.names(Data) %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


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
Aobj<-contMap(pruned.tree,neoplasia,plot=FALSE,res=200, type = "fan")

## invert color scale so red is highest trait value and blue is lowest
Aobj<-setMap(Aobj,invert=TRUE)

## plot 
s<-plot(Aobj,
        fsize=c(1,0.8),leg.txt="neoplasia prevalence",lwd = 2.5,margin = 20)
title(main="C", adj = .05, line = -1)

## plot multiple phylos on to one image

layout(matrix(c(1,1,1,2,2,3), 6,1, byrow = TRUE))
plot(Sobj,
     fsize=c(1,0.8),leg.txt="neoplasia prevalence")
title(main="A", cex.main=3, line = -2, adj=.05)
plot(Mobj,
    fsize=c(1,0.8),leg.txt="neoplasia prevalence")
title(main="B", cex.main=3, line = -2, adj=.05)
plot(Aobj,
    fsize=c(1,0.8),leg.txt="neoplasia prevalence")
title(main="C", cex.main=3, line = -2, adj=.05)




