library(tidyverse)
library(phytools)
library(OUwie)
library(ape)
library(ggpubr)
library(gridGraphics)
library(patchwork)
library(cowplot)
library(devtools)
library(ggtree)
library(plotrix)


## read data and tree frosm file
##Mammals


ogData<-read.csv(file="min20-2022.05.16.csv")
Data <- ogData[,c(4,9),drop=FALSE] 
orders<-ogData[,c(6,9),drop=FALSE] 
Data[Data < 0] <-NA
Data <- na.omit(Data)
fish_species <- data.frame(
  Class = c("Agnatha", "Chondrichthyes", "Osteichthyes"),
  Species = c("Petromyzon marinus", "Carcharodon carcharias", "Gadus morhua")
)
Data<-rbind(Data,fish_species)



tree <- read.tree("superTree.nwk")


#tree<-read.tree(file="subtree.nwk")


Data$Species <- gsub(" ", "_", Data$Species)

orders$Species <- gsub(" ", "_", orders$Species)


pruneTipsWithNumbers <- function(tree) {
  # Identify tips with numbers in their labels using regular expression
  tips_with_numbers <- grep("[0-9]", tree$tip.label)
  
  # If there are any tips with numbers, prune them from the tree
  if (length(tips_with_numbers) > 0) {
    tree <- drop.tip(tree, tree$tip.label[tips_with_numbers])
  }
  
  return(tree)
}

pruneSpecificTips <- function(tree) {
  # Identify tips with numbers, "_sp", or "_cf" in their labels using regular expression
  tips_to_prune <- grep("[0-9]|_sp|_cf", tree$tip.label)
  
  # If there are any tips to prune, remove them from the tree
  if (length(tips_to_prune) > 0) {
    tree <- drop.tip(tree, tree$tip.label[tips_to_prune])
  }
  
  return(tree)
}

pruned.tree <- pruneTipsWithNumbers(tree)
pruned.tree <- pruneSpecificTips(pruned.tree)

name.check(pruned.tree,Data$Species)


pruned.tree$tip.label <- gsub("'", "",pruned.tree$tip.label)
pruned.tree$tip.label

if("Chelodina_mccordi" %in% pruned.tree$tip.label) {
  print("Carcharodon carcharias is in the tree's tip labels.")
} else {
  print("Carcharodon carcharias is not in the tree's tip labels.")
}


#plot(tree, ftype = "off")


inBigtree<-drop.tip(
  pruned.tree, setdiff(
    pruned.tree$tip.label, Data$Species))


speciesMatch<-inBigtree$tip.label


commonAnc<-findMRCA(pruned.tree,speciesMatch)
descendants<-getDescendants(pruned.tree,commonAnc)
pruned.tree<-keep.tip(pruned.tree, descendants)


Data <- Data[Data$Species %in% speciesMatch, ]
orderData<-left_join(Data, orders, by = "Species")



speciesInTree <- pruned.tree$tip.label

# Species in 'Data' that do not match the species in the tree
nonMatchingSpecies <- setdiff(speciesInTree, Data$Species)

# You can print or use 'nonMatchingSpecies' as needed
nonMatchingSpecies

#node_of_interest<-findMRCA(tree, Data$Species)

#pruned.tree <- extract.clade(tree, node = node_of_interest)

#write.tree(pruned.tree, file = "pruned.tree.nwk")

Data <- Data %>%
  filter(!(Class %in% c("Agnatha", "Chondrichthyes", "Osteichthyes")))


mammals<-filter(Data, Class == "Mammalia")
aves<-filter(Data, Class == "Aves")
reptilia<-filter(Data, Class == "Reptilia")
amph<-filter(Data, Class == "Amphibia")

nonMatchingDataFrame <- data.frame(
  Species = nonMatchingSpecies,
  Class = rep("noClade", length(nonMatchingSpecies))
)


combinedDataFrame <- rbind(mammals, aves, reptilia, amph, nonMatchingDataFrame)

classes<-unique(combinedDataFrame$Class)

colors<-c("transparent","#3B4992FF","#EE0000FF","#008280FF","gold","#A20056FF","#5F559BFF")

names(colors)<-0:6

#color code tree branches

# pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, mammals$Species),state="1", anc = "0")
# pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, aves$Species),state="2")
# pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, reptilia$Species),state="3")
# pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, amph$Species),state="4")

getTipIndices <- function(tips, tree) {
  nodes <- sapply(tips, function(x) which(tree$tip.label == x))
  return(nodes)
}

#branch<-getTipIndices("'Epipedobates_machalilla'", pruned.tree)
#paintBranches(pruned.tree, edge=branch, state=black_color_code)
# Assuming 'colors' vector is already defined and includes a black color with code "0"
paintTips <- function(pruned.tree, data_species) {
  # Extract all tip labels from the pruned.tree
  all_tips <- pruned.tree$tip.label
  
  # Loop through each tip label
  for (tip in all_tips) {
    # Check if the tip is not in data_species
    if ((tip %in% data_species)) {
      print(paste("Painting tip:", tip))
      # Get the index of the tip in the tree
      tip_index <- which(pruned.tree$tip.label == tip)
      if (tip %in% mammals$Species){
      # Paint the corresponding branch
        pruned.tree <- paintBranches(pruned.tree, edge=tip_index, state="mammals", anc.state = "0")}
      else if(tip %in% aves$Species){
        pruned.tree <- paintBranches(pruned.tree, edge=tip_index, state="aves", anc.state = "0")}
      else if(tip %in% reptilia$Species){
        pruned.tree <- paintBranches(pruned.tree, edge=tip_index, state="reptilia", anc.state = "0")}
      else if(tip %in% amph$Species){
        pruned.tree <- paintBranches(pruned.tree, edge=tip_index, state="amph", anc.state = "0")}
    }
  }
  
  return(pruned.tree)
}


# Apply the function to paint tips black if they are not in Data$Species
pruned.tree <- paintTips(pruned.tree, Data$Species)

pruned.tree<-reorder.phylo(pruned.tree, order = "cladewise")

#write.tree(pruned.tree, file = "painttree.nwk")

#pruned.tree<-read.tree("painttree.nwk")

# par(fg="black")
# plotSimmap(pruned.tree,colors,type="fan",lwd=0.5,ftype="off", part = 0.5)
# paint<-read.tree(file = "treepaint.tre")

colors <- c("grey87", "transparent", "transparent", "transparent", "transparent", "transparent", "transparent")

# Assign names to the colors based on their index
names(colors)<-c(0, "mammals",'aves','reptilia',"amph")
par(fg="transparent")
plot(pruned.tree,colors,ftype="off",fsize=0.8,lwd=1,type = "fan", part = 0.5)




colors<-c("transparent","#3B4992FF","#EE0000FF","#008280FF","gold","#A20056FF","#5F559BFF")

names(colors)<-c(0, "mammals",'aves','reptilia',"amph")
par(lty="solid",fg="black")

plot(pruned.tree,colors,ftype="off",fsize=0.4,lwd=1, type = "fan",add = TRUE,part = 0.5)
par(lty="solid",fg="black")
legend("topleft",c("Mammals","Aves", "Reptiles", "Amphibians"),
       lwd=1,col=colors[2:5],
       bty="n")
  




rodentia<-filter(orderData, is.element(Orders, c("Rodentia")))
primates<-filter(orderData, is.element(Orders, c("Primates")))
chiroptera<-filter(orderData, is.element(Orders, c("Chiroptera")))
cetacea<-filter(orderData, is.element(Orders, c("Cetacea")))
artio<-filter(orderData, is.element(Orders, c("Artiodactyla")))
cetartio<-rbind(cetacea,artio)
carnivora<-filter(orderData, is.element(Orders, c("Carnivora")))
diprotodontia<-filter(orderData, is.element(Orders, c("Diprotodontia")))
aves<-filter(Data, Class == "Aves")
squamata<-filter(orderData, is.element(Orders, c("Squamata")))
anura<-filter(orderData, is.element(Orders, c("Anura")))
caudata<-filter(orderData, is.element(Orders, c("Caudata")))

species_list <- list(
  Agnatha = c("Petromyzon_marinus", "Eudontomyzon_mariae", "Ichthyomyzon_castaneus"),
  Chondrichthyes = c("Carcharodon_carcharias", "Sphyrna_lewini", "Mobula_mobular"),
  Osteichthyes = c("Salmo_salar", "Gadus_morhua", "Anguilla_anguilla")
)
species_list$Agnatha %in% pruned.tree$tip.label
species_list$Chondrichthyes %in% pruned.tree$tip.label
species_list$Osteichthyes %in% pruned.tree$tip.label


findMRCA(pruned.tree, species_list$Agnatha)


par(fg="#222222")
arc.cladelabels(text="Rodentia",cex = 1.05,node=findMRCA(pruned.tree, rodentia$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)

par(fg="#222222")
arc.cladelabels(text="Cetartiodactyla",cex = 1.05,node=findMRCA(pruned.tree, cetartio$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)
par(fg="#222222")
arc.cladelabels(text="Carnivora",cex = 1.05,node=findMRCA(pruned.tree, carnivora$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)


par(fg="#222222")
arc.cladelabels(text=" Chiroptera",cex = 1.05,node=findMRCA(pruned.tree, chiroptera$Species),ln.offset=1.05,lab.offset=1.1,mark.node =FALSE)

par(fg="#222222")
arc.cladelabels(text="Primates",cex = 1.05,node=findMRCA(pruned.tree, primates$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)

par(fg="#222222")
arc.cladelabels(text="Diprotodontia",cex = 1.05,node=findMRCA(pruned.tree, diprotodontia$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)


par(fg="#222222")
arc.cladelabels(text="Aves",cex = 1.05,node=findMRCA(pruned.tree, aves$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE)


par(fg="#222222")
arc.cladelabels(text="Squamata",cex = 1.05,node=findMRCA(pruned.tree, squamata$Species),ln.offset=1.05,lab.offset=1.1,mark.node=FALSE,mark.node=FALSE)

par(fg="#222222")
arc.cladelabels(text="Anura",cex = 1.05,node=findMRCA(pruned.tree, anura$Species),ln.offset=1.05,lab.offset=1.1, mark.node=FALSE)
par(fg="#222222")
arc.cladelabels(text="Caudata",cex = 1.05,node=findMRCA(pruned.tree, caudata$Species),ln.offset=1.05,lab.offset=1.1, mark.node=FALSE)
par(fg="#222222")
arc.cladelabels(text="Agnatha",cex = 1.05,node=findMRCA(pruned.tree, species_list$Agnatha),ln.offset=1.05,lab.offset=1.1, mark.node=FALSE)



par(fg="#222222")
arc.cladelabels(text="Chondrichthyes",cex = 1.05,node=findMRCA(pruned.tree, species_list$Chondrichthyes),ln.offset=1.05,lab.offset=1.1, mark.node=FALSE)

par(fg="#222222")
arc.cladelabels(text="Osteichthyes",cex = 1.05,node=findMRCA(pruned.tree, species_list$Osteichthyes),ln.offset=1.05,lab.offset=1.1, mark.node=FALSE)
