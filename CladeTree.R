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
#cladelabels(text="Eulipotyphla",node=which(pruned.tree$tip.label=="Four-toed_hedgehog"))
#cladelabels(text="Hyracoidea",cex=.8,node=which(pruned.tree$tip.label=="Rock_hyrax"))
#cladelabels(text="Lagomorpha",cex=.8,node=which(pruned.tree$tip.label=="Domestic_rabbit"))
#cladelabels(text="Perissodactyla",cex = .8,node=which(pruned.tree$tip.label=="Grevys_zebra"))
cladelabels(text="Primates",cex = .8,node=findMRCA(pruned.tree, primates$Species))
#cladelabels(text="Proboscidea",cex = .8,node=which(pruned.tree$tip.label=="Asian_elephant"))
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






#read for Sauropsids
Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
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
fit

## confidence intervals
fit$CI[1,]

## projection of the reconstruction onto the edges of the tree
Sobj<-contMap(pruned.tree,neoplasia,plot=FALSE,res=200, type = "fan")

## invert color scale so red is highest trait value and blue is lowest
Sobj<-setMap(Sobj,invert=TRUE)

## plot 
s<-plot(Sobj,
           fsize=c(1,0.8),leg.txt="neoplasia prevalence",lwd = 2.5, ftype="off")
title(main="B", adj = .05, line = -1)

#order subset
accipitriformes<-filter(specData, is.element(Orders, c("Accipitriformes")))
anseriformes<-filter(specData, is.element(Orders, c("Anseriformes")))
bucerotiformes<-filter(specData, is.element(Orders, c("Bucerotiformes")))
charadriiformes<-filter(specData, is.element(Orders, c("Charadriiformes")))
caprimulgiformes<-filter(specData, is.element(Orders, c("Caprimulgiformes")))
ciconiiformes<-filter(specData, is.element(Orders, c("Ciconiiformes")))
columbiformes<-filter(specData, is.element(Orders, c("Columbiformes")))
crocodilia<-filter(specData, is.element(Orders, c("Crocodilia")))
galliformes<-filter(specData, is.element(Orders, c("Galliformes")))
passeriformes<-filter(specData, is.element(Orders, c("Passeriformes")))
pelecaniformes<-filter(specData, is.element(Orders, c("Pelecaniformes")))
piciformes<-filter(specData, is.element(Orders, c("Piciformes")))
psittaciformes<-filter(specData, is.element(Orders, c("Psittaciformes")))
rheiformes<-filter(specData, is.element(Orders, c("Rheiformes")))
sphenisciformes<-filter(specData, is.element(Orders, c("Sphenisciformes")))
squamata<-filter(specData, is.element(Orders, c("Squamata")))
strigiformes<-filter(specData, is.element(Orders, c("Strigiformes")))
testudines<-filter(specData, is.element(Orders, c("Testudines")))

#Family subset
Accipitridae<-filter(specData, is.element(Family, c("Accipitridae")))
Agamidae<-filter(specData, is.element(Family, c("Agamidae")))
Alligatoridae<-filter(specData, is.element(Family, c("Alligatoridae")))
Anatidae<-filter(specData, is.element(Family, c("Anatidae")))
Anguidae<-filter(specData, is.element(Family, c("Anguidae")))
Ardeidae<-filter(specData, is.element(Family, c("Ardeidae")))
Boidae<-filter(specData, is.element(Family, c("Boidae")))
Bucerotidae<-filter(specData, is.element(Family, c("Bucerotidae")))
Cacatuidae<-filter(specData, is.element(Family, c("Cacatuidae")))
Chamaeleonidae<-filter(specData, is.element(Family, c("Chamaeleonidae")))
Ciconiidae<-filter(specData, is.element(Family, c("Ciconiidae")))
Colubridae<-filter(specData, is.element(Family, c("Colubridae")))
Columbidae<-filter(specData, is.element(Family, c("Columbidae")))
Corytophanidae<-filter(specData, is.element(Family, c("Corytophanidae")))
Crotaphytidae<-filter(specData, is.element(Family, c("Crotaphytidae")))
Dactyloidae<-filter(specData, is.element(Family, c("Dactyloidae")))
Elapidae<-filter(specData, is.element(Family, c("Elapidae")))
Emberizidae<-filter(specData, is.element(Family, c("Emberizidae")))
Estrildidae<-filter(specData, is.element(Family, c("Estrildidae")))
Eurypygidae<-filter(specData, is.element(Family, c("Eurypygidae")))
Falconidae<-filter(specData, is.element(Family, c("Falconidae")))
Fringillidae<-filter(specData, is.element(Family, c("Fringillidae")))
Gekkonidae<-filter(specData, is.element(Family, c("Gekkonidae")))
Helodermatidae<-filter(specData, is.element(Family, c("Helodermatidae")))
Homalopsidae<-filter(specData, is.element(Family, c("Homalopsidae")))
Iguanidae<-filter(specData, is.element(Family, c("Iguanidae")))
Laridae<-filter(specData, is.element(Family, c("Laridae")))
Leiothrichidae<-filter(specData, is.element(Family, c("Leiothrichidae")))
Loriidae<-filter(specData, is.element(Family, c("Loriidae")))
Momotidae<-filter(specData, is.element(Family, c("Momotidae")))
Muscicapidae<-filter(specData, is.element(Family, c("Muscicapidae")))
Musophagidae<-filter(specData, is.element(Family, c("Musophagidae")))
Numididae<-filter(specData, is.element(Family, c("Numididae")))
Odontophoridae<-filter(specData, is.element(Family, c("Odontophoridae")))
Pelecanidae<-filter(specData, is.element(Family, c("Pelecanidae")))
Phasianidae<-filter(specData, is.element(Family, c("Phasianidae")))
Phrynosomatidae<-filter(specData, is.element(Family, c("Phrynosomatidae")))
Pittidae<-filter(specData, is.element(Family, c("Pittidae")))
Ploceidae<-filter(specData, is.element(Family, c("Ploceidae")))
Podargidae<-filter(specData, is.element(Family, c("Podargidae")))
Psittacidae<-filter(specData, is.element(Family, c("Psittacidae")))
Psittaculidae<-filter(specData, is.element(Family, c("Psittaculidae")))
Pythonidae<-filter(specData, is.element(Family, c("Pythonidae")))
Rallidae<-filter(specData, is.element(Family, c("Rallidae")))
Ramphastidae<-filter(specData, is.element(Family, c("Ramphastidae")))
Rheidae<-filter(specData, is.element(Family, c("Rheidae")))
Scincidae<-filter(specData, is.element(Family, c("Scincidae")))
Scopidae<-filter(specData, is.element(Family, c("Scopidae")))
Spheniscidae<-filter(specData, is.element(Family, c("Spheniscidae")))
Strigidae<-filter(specData, is.element(Family, c("Strigidae")))
Sturnidae<-filter(specData, is.element(Family, c("Sturnidae")))
Testudinidae<-filter(specData, is.element(Family, c("Testudinidae")))
Thraupidae<-filter(specData, is.element(Family, c("Thraupidae")))
Threskiornithidae<-filter(specData, is.element(Family, c("Threskiornithidae")))
Tinamidae<-filter(specData, is.element(Family, c("Tinamidae")))
Tropiduridae<-filter(specData, is.element(Family, c("Tropiduridae")))
Turdidae<-filter(specData, is.element(Family, c("Turdidae")))
Varanidae<-filter(specData, is.element(Family, c("Varanidae")))





#order lines
par(fg="#000000")
cladelabels(text="Accipitriformes",cex = .8,node=findMRCA(pruned.tree, accipitriformes$Species ))
cladelabels(text="Anseriformes",cex = .8,node=findMRCA(pruned.tree, anseriformes$Species ) )
cladelabels(text="Bucerotiformes",cex = .8,node=findMRCA(pruned.tree, bucerotiformes$Species ))
cladelabels(text="Charadriiformes",cex = .8,node=findMRCA(pruned.tree, charadriiformes$Species ))
cladelabels(text="Caprimulgiformes",cex = .8,node=findMRCA(pruned.tree, caprimulgiformes$Species ) )
cladelabels(text="Ciconiiformes",cex = .8,node=findMRCA(pruned.tree, ciconiiformes$Species ) )
cladelabels(text="Columbiformes",cex = .8,node=findMRCA(pruned.tree, columbiformes$Species ) )
cladelabels(text="Crocodilia",cex = .8,node=findMRCA(pruned.tree, crocodilia$Species ) )
cladelabels(text="Galliformes",cex = .8,node=findMRCA(pruned.tree, galliformes$Species ) )
cladelabels(text="Pelecaniformes",cex = .8,node=findMRCA(pruned.tree, pelecaniformes$Species ) )
cladelabels(text="Passeriformes",cex = .8,node=findMRCA(pruned.tree, passeriformes$Species ))
cladelabels(text="piciformes",cex = .8,node=findMRCA(pruned.tree, piciformes$Species ) )
cladelabels(text="Psittaciformes",cex = .8,node=findMRCA(pruned.tree, psittaciformes$Species ))
cladelabels(text="Rheiformes",cex = .8,node=findMRCA(pruned.tree, rheiformes$Species ))
cladelabels(text="Sphenisciformes",cex = .8,node=findMRCA(pruned.tree, sphenisciformes$Species ))
cladelabels(text="Squamata",cex = .8,node=findMRCA(pruned.tree, squamata$Species ))
cladelabels(text="Strigiformes",cex = .8,node=findMRCA(pruned.tree, strigiformes$Species ))
cladelabels(text="Testudines",cex = .8,node=findMRCA(pruned.tree, testudines$Species ))



#family lines
par(fg="#808080")
cladelabels(text="Accipitridae",cex = .8,node=findMRCA(pruned.tree, Accipitridae$Species))
cladelabels(text="Agamidae",cex = .8,node=findMRCA(pruned.tree, Agamidae$Species))
cladelabels(text="Alligatoridae",cex = .8,node=findMRCA(pruned.tree, Alligatoridae$Species))
cladelabels(text="Anatidae",cex = .8,node=findMRCA(pruned.tree, Anatidae$Species))
cladelabels(text="Anguidae",cex = .8,node=findMRCA(pruned.tree, Anguidae$Species))
cladelabels(text="Ardeidae",cex = .8,node=findMRCA(pruned.tree, Ardeidae$Species))
cladelabels(text="Boidae",cex = .8,node=findMRCA(pruned.tree, Boidae$Species))
cladelabels(text="Bucerotidae",cex = .8,node=findMRCA(pruned.tree, Bucerotidae$Species))
cladelabels(text="Cacatuidae",cex = .8,node=findMRCA(pruned.tree, Cacatuidae$Species))
cladelabels(text="Chamaeleonidae",cex = .8,node=findMRCA(pruned.tree, Chamaeleonidae$Species))
cladelabels(text="Ciconiidae",cex = .8,node=findMRCA(pruned.tree, Ciconiidae$Species))
cladelabels(text="Colubridae",cex = .8,node=findMRCA(pruned.tree, Colubridae$Species))
cladelabels(text="Columbidae",cex = .8,node=findMRCA(pruned.tree, Columbidae$Species))
cladelabels(text="Corytophanidae",cex = .8,node=findMRCA(pruned.tree, Corytophanidae$Species))
cladelabels(text="Crotaphytidae",cex = .8,node=findMRCA(pruned.tree, Crotaphytidae$Species))
cladelabels(text="Dactyloidae",cex = .8,node=findMRCA(pruned.tree, Dactyloidae$Species))
cladelabels(text="Elapidae",cex = .8,node=findMRCA(pruned.tree, Elapidae$Species))
cladelabels(text="Emberizidae",cex = .8,node=findMRCA(pruned.tree, Emberizidae$Species))
cladelabels(text="Estrildidae",cex = .8,node=findMRCA(pruned.tree, Estrildidae$Species))
cladelabels(text="Eurypygidae",cex = .8,node=findMRCA(pruned.tree, Eurypygidae$Species))
cladelabels(text="Falconidae",cex = .8,node=findMRCA(pruned.tree, Falconidae$Species))
cladelabels(text="Fringillidae",cex = .8,node=findMRCA(pruned.tree, Fringillidae$Species))
cladelabels(text="Gekkonidae",cex = .8,node=findMRCA(pruned.tree, Gekkonidae$Species))
cladelabels(text="Helodermatidae",cex = .8,node=findMRCA(pruned.tree, Helodermatidae$Species))
cladelabels(text="Homalopsidae",cex = .8,node=findMRCA(pruned.tree, Homalopsidae$Species))
cladelabels(text="Iguanidae",cex = .8,node=findMRCA(pruned.tree, Iguanidae$Species))
cladelabels(text="Laridae",cex = .8,node=findMRCA(pruned.tree, Laridae$Species))
cladelabels(text="Leiothrichidae",cex = .8,node=findMRCA(pruned.tree, Leiothrichidae$Species))
cladelabels(text="Loriidae",cex = .8,node=findMRCA(pruned.tree, Loriidae$Species))
cladelabels(text="Momotidae",cex = .8,node=findMRCA(pruned.tree, Momotidae$Species))
cladelabels(text="Muscicapidae",cex = .8,node=findMRCA(pruned.tree, Muscicapidae$Species))
cladelabels(text="Musophagidae",cex = .8,node=findMRCA(pruned.tree, Musophagidae$Species))
cladelabels(text="Numididae",cex = .8,node=findMRCA(pruned.tree, Numididae$Species))
cladelabels(text="Odontophoridae",cex = .8,node=findMRCA(pruned.tree, Odontophoridae$Species))
cladelabels(text="Pelecanidae",cex = .8,node=findMRCA(pruned.tree, Pelecanidae$Species))
cladelabels(text="Phasianidae",cex = .8,node=findMRCA(pruned.tree, Phasianidae$Species))
cladelabels(text="Phrynosomatidae",cex = .8,node=findMRCA(pruned.tree, Phrynosomatidae$Species))
cladelabels(text="Pittidae",cex = .8,node=findMRCA(pruned.tree, Pittidae$Species))
cladelabels(text="Ploceidae",cex = .8,node=findMRCA(pruned.tree, Ploceidae$Species))
cladelabels(text="Podargidae",cex = .8,node=findMRCA(pruned.tree, Podargidae$Species))
cladelabels(text="Psittacidae",cex = .8,node=findMRCA(pruned.tree, Psittacidae$Species))
cladelabels(text="Psittaculidae",cex = .8,node=findMRCA(pruned.tree, Psittaculidae$Species))
cladelabels(text="Pythonidae",cex = .8,node=findMRCA(pruned.tree, Pythonidae$Species))
cladelabels(text="Rallidae",cex = .8,node=findMRCA(pruned.tree, Rallidae$Species))
cladelabels(text="Ramphastidae",cex = .8,node=findMRCA(pruned.tree, Ramphastidae$Species))
cladelabels(text="Rheidae",cex = .8,node=findMRCA(pruned.tree, Rheidae$Species))
cladelabels(text="Scincidae",cex = .8,node=findMRCA(pruned.tree, Scincidae$Species))
cladelabels(text="Scopidae",cex = .8,node=findMRCA(pruned.tree, Scopidae$Species))
cladelabels(text="Spheniscidae",cex = .8,node=findMRCA(pruned.tree, Spheniscidae$Species))
cladelabels(text="Strigidae",cex = .8,node=findMRCA(pruned.tree, Strigidae$Species))
cladelabels(text="Sturnidae",cex = .8,node=findMRCA(pruned.tree, Sturnidae$Species))
cladelabels(text="Testudinidae",cex = .8,node=findMRCA(pruned.tree, Testudinidae$Species))
cladelabels(text="Thraupidae",cex = .8,node=findMRCA(pruned.tree, Thraupidae$Species))
cladelabels(text="Threskiornithidae",cex = .8,node=findMRCA(pruned.tree, Threskiornithidae$Species))
cladelabels(text="Tinamidae",cex = .8,node=findMRCA(pruned.tree, Tinamidae$Species))
cladelabels(text="Tropiduridae",cex = .8,node=findMRCA(pruned.tree, Tropiduridae$Species))
cladelabels(text="Turdidae",cex = .8,node=findMRCA(pruned.tree, Turdidae$Species))
cladelabels(text="Varanidae",cex = .8,node=findMRCA(pruned.tree, Varanidae$Species))


#read for Amphibians
Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Amphibia")))
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
        fsize=c(1,0.8),leg.txt="neoplasia prevalence",lwd = 2.5,margin = 20, ftype="off")
title(main="C", adj = .05, line = -1)

#subset orders
anura<-filter(specData, is.element(Orders, c("Anura")))
caudata<-filter(specData, is.element(Orders, c("Caudata")))
gymnophiona<-filter(specData, is.element(Orders, c("Gymnophiona")))



#label orders
par(fg="#000000")
cladelabels(text="Anura",cex = .8,node=findMRCA(pruned.tree, anura$Species))
cladelabels(text="Caudata",cex = .8,node=findMRCA(pruned.tree, caudata$Species))
cladelabels(text="Gymnophiona",cex=.8,node=which(pruned.tree$tip.label=="Gaboon_caecilian"))



#Family subset
Alytidae<-filter(specData, is.element(Family, c("Alytidae")))
Ambystomatidae<-filter(specData, is.element(Family, c("Ambystomatidae")))
Aromobatidae<-filter(specData, is.element(Family, c("Aromobatidae")))
Bufonidae<-filter(specData, is.element(Family, c("Bufonidae")))
Caeciliidae<-filter(specData, is.element(Family, c("Caeciliidae")))
Ceratobatrachidae<-filter(specData, is.element(Family, c("Ceratobatrachidae")))
Ceratophryidae<-filter(specData, is.element(Family, c("Ceratophryidae")))
Dendrobatidae<-filter(specData, is.element(Family, c("Dendrobatidae")))
Dermophiidae<-filter(specData, is.element(Family, c("Dermophiidae")))
Hylidae<-filter(specData, is.element(Family, c("Hylidae")))
Leiopelmatidae<-filter(specData, is.element(Family, c("Leiopelmatidae")))
Leptodactylidae<-filter(specData, is.element(Family, c("Leptodactylidae")))
Mantellidae<-filter(specData, is.element(Family, c("Mantellidae")))
Megophryidae<-filter(specData, is.element(Family, c("Megophryidae")))
Microhylidae<-filter(specData, is.element(Family, c("Microhylidae")))
Pipidae<-filter(specData, is.element(Family, c("Pipidae")))
Ranidae<-filter(specData, is.element(Family, c("Ranidae")))
Rhacophoridae<-filter(specData, is.element(Family, c("Rhacophoridae")))
Salamandridae<-filter(specData, is.element(Family, c("Salamandridae")))


#family lines
par(fg="#808080")
cladelabels(text="Alytidae",cex = .8,node=findMRCA(pruned.tree, Alytidae$Species))
cladelabels(text="Ambystomatidae",cex = .8,node=findMRCA(pruned.tree, Ambystomatidae$Species))
cladelabels(text="Aromobatidae",cex = .8,node=findMRCA(pruned.tree, Aromobatidae$Species))
cladelabels(text="Bufonidae",cex = .8,node=findMRCA(pruned.tree, Bufonidae$Species))
cladelabels(text="Caeciliidae",cex = .8,node=findMRCA(pruned.tree, Caeciliidae$Species))
cladelabels(text="Ceratobatrachidae",cex = .8,node=findMRCA(pruned.tree, Ceratobatrachidae$Species))
cladelabels(text="Ceratophryidae",cex = .8,node=findMRCA(pruned.tree, Ceratophryidae$Species))
cladelabels(text="Dendrobatidae",cex = .8,node=findMRCA(pruned.tree, Dendrobatidae$Species))
cladelabels(text="Dermophiidae",cex = .8,node=findMRCA(pruned.tree, Dermophiidae$Species))
cladelabels(text="Hylidae",cex = .8,node=findMRCA(pruned.tree, Hylidae$Species))
cladelabels(text="Leiopelmatidae",cex = .8,node=findMRCA(pruned.tree, Leiopelmatidae$Species))
cladelabels(text="Leptodactylidae",cex = .8,node=findMRCA(pruned.tree, Leptodactylidae$Species))
cladelabels(text="Mantellidae",cex = .8,node=findMRCA(pruned.tree, Mantellidae$Species))
cladelabels(text="Megophryidae",cex = .8,node=findMRCA(pruned.tree, Megophryidae$Species))
cladelabels(text="Microhylidae",cex = .8,node=findMRCA(pruned.tree, Microhylidae$Species))
cladelabels(text="Pipidae",cex = .8,node=findMRCA(pruned.tree, Pipidae$Species))
cladelabels(text="Ranidae",cex = .8,node=findMRCA(pruned.tree, Ranidae$Species))
cladelabels(text="Rhacophoridae",cex = .8,node=findMRCA(pruned.tree, Rhacophoridae$Species))
cladelabels(text="Salamandridae",cex = .8,node=findMRCA(pruned.tree, Salamandridae$Species))



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




