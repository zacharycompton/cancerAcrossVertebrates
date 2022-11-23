cutData<- filter(specData, is.element(Family, c("	Psittaciformes")))
tree<-read.tree(file="min20Fixed516.nwk")
cutData$Species <- gsub(" ", "_", cutData$Species) 
cutData$common_name<-gsub("_", "", cutData$common_name)
includedSpecies <- cutData$Species

pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
cutData$Keep <- cutData$Species %in% pruned.tree$tip.label
cutData <- cutData[!(cutData$Keep==FALSE),]
rownames(cutData)<-cutData$Species