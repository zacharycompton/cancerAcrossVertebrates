#open libraries
library(tidyverse)
library(phytools)
library(patchwork)
library(tidyverse)
library(cowplot)



#par(mfrow=c(3,1))
#read tree
Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
Data <- Data[,c(6,7,9,10,17,13),drop=FALSE] 
tree<-read.tree(file="min20Fixed516.nwk")
Data[Data < 0] <-NA
Data <- na.omit(Data)

Data$Species <- gsub(" ", "_", Data$Species)
includedSpecies <- Data$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$Species %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]

specData<-Data

rownames(Data)<-Data$Species
Data <- Data[,c(5,6),drop=FALSE] 
#name.check(pruned.tree,Data)


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
     fsize=c(1,0.8),leg.txt="neoplasia prevalence",lwd= 2.5,ftype="off")
title(main="A", adj = .05, line = -1)

specData$Species <- gsub(" ", "_", specData$Species)

#order subsets
rodentia<-filter(specData, is.element(Orders, c("Rodentia")))
artio<-filter(specData, is.element(Orders, c("Artiodactyla")))
carnivora<-filter(specData, is.element(Orders, c("Carnivora")))
cetacea<-filter(specData, is.element(Orders, c("Cetacea")))
chiroptera<-filter(specData, is.element(Orders, c("Chiroptera")))
cingulata<-filter(specData, is.element(Orders, c("Cingulata")))
didelphimorphia<-filter(specData, is.element(Orders, c("Didelphimorphia")))
eulipotyphla<-filter(specData, is.element(Orders, c("Eulipotyphla")))
hyracoidea<-filter(specData, is.element(Orders, c("Hyracoidea")))
lagomorpha<-filter(specData, is.element(Orders, c("Lagomorpha")))
perissodactyla<-filter(specData, is.element(Orders, c("Perissodactyla")))
primates<-filter(specData, is.element(Orders, c("Primates")))
proboscidea<-filter(specData, is.element(Orders, c("Proboscidea")))
diprotodontia<-filter(specData, is.element(Orders, c("Diprotodontia")))

#family subsets
Acrobatidae<-filter(specData, is.element(Family, c("Acrobatidae")))
Ailuridae<-filter(specData, is.element(Family, c("Ailuridae")))
Bathyergidae<-filter(specData, is.element(Family, c("Bathyergidae")))
Bovidae<-filter(specData, is.element(Family, c("Bovidae")))
Callitrichidae<-filter(specData, is.element(Family, c("Callitrichidae")))
Camelidae<-filter(specData, is.element(Family, c("Camelidae")))
Canidae<-filter(specData, is.element(Family, c("Canidae")))
Caviidae<-filter(specData, is.element(Family, c("Caviidae")))
Cebidae<-filter(specData, is.element(Family, c("Cebidae")))
Cercopithecidae<-filter(specData, is.element(Family, c("Cercopithecidae")))
Cervidae<-filter(specData, is.element(Family, c("Cervidae")))
Chinchillidae<-filter(specData, is.element(Family, c("Chinchillidae")))
Cricetidae<-filter(specData, is.element(Family, c("Cricetidae")))
Dasypodidae<-filter(specData, is.element(Family, c("Dasypodidae")))
Delphinidae<-filter(specData, is.element(Family, c("Delphinidae")))
Elephantidae<-filter(specData, is.element(Family, c("Elephantidae")))
Elephantidae<-filter(specData, is.element(Family, c("Elephantidae")))
Equidae<-filter(specData, is.element(Family, c("Equidae")))
Erinaceidae<-filter(specData, is.element(Family, c("Erinaceidae")))
Felidae<-filter(specData, is.element(Family, c("Felidae")))
Giraffidae<-filter(specData, is.element(Family, c("Giraffidae")))
Gliridae<-filter(specData, is.element(Family, c("Gliridae")))
Herpestidae<-filter(specData, is.element(Family, c("Herpestidae")))
Hominidae<-filter(specData, is.element(Family, c("Hominidae")))
Hydrochaeridae<-filter(specData, is.element(Family, c("Hydrochaeridae")))
Lemuridae<-filter(specData, is.element(Family, c("Lemuridae")))
Leporidae<-filter(specData, is.element(Family, c("Leporidae")))
Macropodidae<-filter(specData, is.element(Family, c("Macropodidae")))
Macroscelididae<-filter(specData, is.element(Family, c("Macroscelididae")))
Megalonychidae<-filter(specData, is.element(Family, c("Megalonychidae")))
Muridae<-filter(specData, is.element(Family, c("Muridae")))
Mustelidae<-filter(specData, is.element(Family, c("Mustelidae")))
Octodontidae<-filter(specData, is.element(Family, c("Octodontidae")))
Otariidae<-filter(specData, is.element(Family, c("Otariidae")))
Petauridae<-filter(specData, is.element(Family, c("Petauridae")))
Phascolarctidae<-filter(specData, is.element(Family, c("Phascolarctidae")))
Phocidae<-filter(specData, is.element(Family, c("Phocidae")))
Phocoenidae<-filter(specData, is.element(Family, c("Phocoenidae")))
Phyllostomidae<-filter(specData, is.element(Family, c("Phyllostomidae")))
Pitheciidae<-filter(specData, is.element(Family, c("Pitheciidae")))
Procaviidae<-filter(specData, is.element(Family, c("Procaviidae")))
Procyonidae<-filter(specData, is.element(Family, c("Procyonidae")))
Pteropodidae<-filter(specData, is.element(Family, c("Pteropodidae")))
Sciuridae<-filter(specData, is.element(Family, c("Sciuridae")))
Suidae<-filter(specData, is.element(Family, c("Suidae")))
Tayassuidae<-filter(specData, is.element(Family, c("Tayassuidae")))
Ursidae<-filter(specData, is.element(Family, c("Ursidae")))





#order lines
par(fg="#000000")
cladelabels(tree = pruned.tree,text="Rodentia",cex = .8,node=findMRCA(pruned.tree, rodentia$Species))
cladelabels(text="Artiodactyla",cex = .8,node=findMRCA(pruned.tree, artio$Species))
cladelabels(text="Carnivora",cex = .8,node=findMRCA(pruned.tree, carnivora$Species))
cladelabels(text="Cetacea",cex = .8,node=findMRCA(pruned.tree, cetacea$Species))
cladelabels(text="Chiroptera",cex = .8,node=findMRCA(pruned.tree, chiroptera$Species))
cladelabels(text="Didelphimorphia",cex=.8,node=findMRCA(pruned.tree, didelphimorphia$Species))
cladelabels(text="Eulipotyphla",node=which(pruned.tree$tip.label=="Four-toed_hedgehog"))
cladelabels(text="Hyracoidea",cex=.8,node=which(pruned.tree$tip.label=="Rock_hyrax"))
cladelabels(text="Lagomorpha",cex=.8,node=which(pruned.tree$tip.label=="Domestic_rabbit"))
cladelabels(text="Perissodactyla",cex = .8,node=which(pruned.tree$tip.label=="Grevys_zebra"))
cladelabels(text="Primates",cex = .8,node=findMRCA(pruned.tree, primates$Species))
cladelabels(text="Proboscidea",cex = .8,node=which(pruned.tree$tip.label=="Asian_elephant"))
cladelabels(text="Diprotodontia",cex = .8,node=findMRCA(pruned.tree, diprotodontia$Species))


#family lines
par(fg="#808080")
cladelabels(text="Acrobatidae",cex = .8,node=findMRCA(pruned.tree, Acrobatidae$Species))
cladelabels(text="Ailuridae",cex = .8,node=findMRCA(pruned.tree, Ailuridae$Species))
cladelabels(text="Bathyergidae",cex = .8,node=findMRCA(pruned.tree, Bathyergidae$Species))
cladelabels(text="Bovidae",cex = .8,node=findMRCA(pruned.tree, Bovidae$Species))
cladelabels(text="Callitrichidae",cex = .8,node=findMRCA(pruned.tree, Callitrichidae$Species))
cladelabels(text="Camelidae",cex = .8,node=findMRCA(pruned.tree, Camelidae$Species))
cladelabels(text="Canidae",cex = .8,node=findMRCA(pruned.tree, Canidae$Species))
cladelabels(text="Caviidae",cex = .8,node=findMRCA(pruned.tree, Caviidae$Species))
cladelabels(text="Cebidae",cex = .8,node=findMRCA(pruned.tree, Cebidae$Species))
cladelabels(text="Cercopithecidae",cex = .8,node=findMRCA(pruned.tree, Cercopithecidae$Species))
cladelabels(text="Cervidae",cex = .8,node=findMRCA(pruned.tree, Cervidae$Species))
cladelabels(text="Chinchillidae",cex = .8,node=findMRCA(pruned.tree, Chinchillidae$Species))
cladelabels(text="Cricetidae",cex = .8,node=findMRCA(pruned.tree, Cricetidae$Species))
cladelabels(text="Dasypodidae",cex = .8,node=findMRCA(pruned.tree, Dasypodidae$Species))
cladelabels(text="Delphinidae",cex = .8,node=findMRCA(pruned.tree, Delphinidae$Species))
cladelabels(text="Elephantidae",cex = .8,node=findMRCA(pruned.tree, Elephantidae$Species))
cladelabels(text="Equidae",cex = .8,node=findMRCA(pruned.tree, Equidae$Species))
cladelabels(text="Erinaceidae",cex = .8,node=findMRCA(pruned.tree, Erinaceidae$Species))
cladelabels(text="Felidae",cex = .8,node=findMRCA(pruned.tree, Felidae$Species))
cladelabels(text="Giraffidae",cex = .8,node=findMRCA(pruned.tree, Giraffidae$Species))
cladelabels(text="Gliridae",cex = .8,node=findMRCA(pruned.tree, Gliridae$Species))
cladelabels(text="Herpestidae",cex = .8,node=findMRCA(pruned.tree, Herpestidae$Species))
cladelabels(text="Hominidae",cex = .8,node=findMRCA(pruned.tree, Hominidae$Species))
cladelabels(text="Hydrochaeridae",cex = .8,node=findMRCA(pruned.tree, Hydrochaeridae$Species))
cladelabels(text="Lemuridae",cex = .8,node=findMRCA(pruned.tree, Lemuridae$Species))
cladelabels(text="Leporidae",cex = .8,node=findMRCA(pruned.tree, Leporidae$Species))
cladelabels(text="Macropodidae",cex = .8,node=findMRCA(pruned.tree, Macropodidae$Species))
cladelabels(text="Macroscelididae",cex = .8,node=findMRCA(pruned.tree, Macroscelididae$Species))
cladelabels(text="Megalonychidae",cex = .8,node=findMRCA(pruned.tree, Megalonychidae$Species))
cladelabels(text="Muridae",cex = .8,node=findMRCA(pruned.tree, Muridae$Species))
cladelabels(text="Mustelidae",cex = .8,node=findMRCA(pruned.tree, Mustelidae$Species))
cladelabels(text="Octodontidae",cex = .8,node=findMRCA(pruned.tree, Octodontidae$Species))
cladelabels(text="Otariidae",cex = .8,node=findMRCA(pruned.tree, Otariidae$Species))
cladelabels(text="Petauridae",cex = .8,node=findMRCA(pruned.tree, Petauridae$Species))
cladelabels(text="Phascolarctidae",cex = .8,node=findMRCA(pruned.tree, Phascolarctidae$Species))
cladelabels(text="Phocidae",cex = .8,node=findMRCA(pruned.tree, Phocidae$Species))
cladelabels(text="Phocoenidae",cex = .8,node=findMRCA(pruned.tree, Phocoenidae$Species))
cladelabels(text="Phyllostomidae",cex = .8,node=findMRCA(pruned.tree, Phyllostomidae$Species))
cladelabels(text="Pitheciidae",cex = .8,node=findMRCA(pruned.tree, Pitheciidae$Species))
cladelabels(text="Procaviidae",cex = .8,node=findMRCA(pruned.tree, Procaviidae$Species))
cladelabels(text="Procyonidae",cex = .8,node=findMRCA(pruned.tree, Procyonidae$Species))
cladelabels(text="Pteropodidae",cex = .8,node=findMRCA(pruned.tree, Pteropodidae$Species))
cladelabels(text="Sciuridae",cex = .8,node=findMRCA(pruned.tree, Sciuridae$Species))
cladelabels(text="Suidae",cex = .8,node=findMRCA(pruned.tree, Suidae$Species))
cladelabels(text="Tayassuidae",cex = .8,node=findMRCA(pruned.tree, Tayassuidae$Species))
cladelabels(text="Ursidae",cex = .8,node=findMRCA(pruned.tree, Ursidae$Species))







## read csv file for Sauropsids

Data<-read.csv("min20516.csv",header = TRUE,row.names=9)  
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
row.names(Data)= gsub(" ", "_", row.names(Data))


## prune the tree to match the data
includedSpecies <- row.names(Data)
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
includedSpecies <- row.names(Data)
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




