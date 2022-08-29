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
Data <- Data[,c(6,9,10,17,13),drop=FALSE] 
tree<-read.tree(file="min20Fixed516.nwk")
Data[Data < 0] <-NA
Data <- na.omit(Data)

Data$Species <- gsub(" ", "_", Data$Species)
includedSpecies <- Data$Species
newtips<-str_remove_all(tree$tip.label,"_ott")
newtips<-str_remove_all(newtips,".ott")
newtips<-str_remove_all(newtips,"-ott")
newtips<-str_remove_all(newtips,"[1234567890]")
newtips<-sub('^([^_]+_[^_]+).*', '\\1', newtips)
tree$tip.label <- newtips
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$Species %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]

specData<-Data

rownames(Data)<-Data$Species
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)


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

par(fg="#3B4992FF")
cladelabels(tree = pruned.tree,text="Rodentia",cex = .8,node=findMRCA(pruned.tree, rodentia$Species))
nodelabels(node=findMRCA(pruned.tree, rodentia$Species), frame = "circle", bg="#3B4992FF", text=" ", )
par(fg="#EE0000FF")
cladelabels(text="Artiodactyla",cex = .8,node=findMRCA(pruned.tree, artio$Species))
nodelabels(node=findMRCA(pruned.tree, artio$Species), frame = "circle", bg="#EE0000FF", text=" ", )
par(fg="#008280FF")
cladelabels(text="Carnivora",cex = .8,node=findMRCA(pruned.tree, carnivora$Species))
nodelabels(node=findMRCA(pruned.tree, carnivora$Species), frame = "circle", bg="#008280FF", text=" ", )
#cladelabels(text="Cetacea",cex = .8,node=findMRCA(pruned.tree, cetacea$Species),ln.offset=1.7,lab.offset=1.75)
#cladelabels(text=" Chiroptera",cex = .8,node=findMRCA(pruned.tree, chiroptera$Species),ln.offset=1.7,lab.offset=1.75)
#cladelabels(text="Didelphimorphia",cex=.5,node=findMRCA(pruned.tree, didelphimorphia$Species),ln.offset=1.7,lab.offset=1.8, orientation = "horizontal")
#cladelabels(text="Eulipotyphla",node=which(pruned.tree$tip.label=="Four-toed_hedgehog"), orientation="horizontal",ln.offset=1.7,lab.offset=1.75)
#cladelabels(text="Hyracoidea",fsize=.4,node=65,orientation="horizontal",ln.offset=1.45,lab.offset=1.5,mark.node=FALSE)
#cladelabels(text="Lagomorpha",cex.sub=.1,node=which(pruned.tree$tip.label=="Domestic_rabbit"),orientation="horizontal",ln.offset=1.45,lab.offset=1.45, mark.node=FALSE)
#cladelabels(text="Perissodactyla",cex = .8,node=which(pruned.tree$tip.label=="Grevys_zebra"),orientation="horizontal",ln.offset=1.45,lab.offset=1.5, mark.node=FALSE)
par(fg="#A20056FF")
cladelabels(text="Primates",cex = .8,node=findMRCA(pruned.tree, primates$Species))
nodelabels(node=findMRCA(pruned.tree, primates$Species), frame = "circle", bg="#A20056FF", text=" ", )
#cladelabels(text="Proboscidea",cex = .8,node=which(pruned.tree$tip.label=="Asian_elephant"),orientation="horizontal",ln.offset=1.45,lab.offset=1.5, mark.node=FALSE)
par(fg="#5F559BFF")
cladelabels(text="Diprotodontia",cex = .8,node=findMRCA(pruned.tree, diprotodontia$Species))
nodelabels(node=findMRCA(pruned.tree, diprotodontia$Species), frame = "circle", bg="#5F559BFF", text=" ", )



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




