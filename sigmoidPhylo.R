###function for the sigmoid phylo
sigmoidPhylogram<-function(tree,...){
  ## b=5,m1=0.01,m2=0.5,v=1
  b<-if(hasArg(b)) list(...)$b else 5
  m1<-if(hasArg(m1)) list(...)$m1 else 0.01
  m2<-if(hasArg(m2)) list(...)$m2 else 0.5
  v<-if(hasArg(v)) list(...)$v else 1
  if(inherits(tree,"simmap")){
    if(hasArg(colors)) colors<-list(...)$colors
    else {
      ss<-sort(unique(c(getStates(tree,"nodes"),
                        getStates(tree,"tips"))))
      mm<-length(ss)
      colors<-setNames(
        colorRampPalette(palette()[1:min(8,mm)])(mm),
        ss)
    }
  } else if(inherits(tree,"phylo")) {
    if(hasArg(color)) colors<-setNames(list(...)$color,"1")
    else colors<-setNames(par()$fg,"1")
    tree<-paintSubTree(tree,Ntip(tree)+1,"1")
  }
  if(hasArg(res)) res<-list(...)$res
  else res<-199
  if(hasArg(use.edge.length))
    use.edge.length<-list(...)$use.edge.length
  else use.edge.length<-TRUE
  if(!use.edge.length){
    if(hasArg(power)) power<-list(...)$power
    else power<-1
    tree<-compute.brlen.simmap(tree,power=power)
  }
  if(hasArg(lwd)) lwd<-list(...)$lwd
  else lwd<-2
  h<-max(nodeHeights(tree))
  args<-list(...)
  args$power<-NULL
  args$res<-NULL
  args$colors<-NULL
  args$b<-NULL
  args$m1<-NULL
  args$m2<-NULL
  args$v<-NULL
  args$tree<-tree
  args$color<-"transparent"
  dev.hold()
  do.call(plotTree,args)
  pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## Yt<-A+(K-A)/((C+exp(-B*(t-M)))^(1/v))
  sigmoid<-function(x,.A=A,.K=K,.C=C,.B=B,.M=M,.v=v)
    return(.A+(.K-.A)/((.C+exp(-.B*(x-.M)))^(1/.v)))
  for(i in 1:nrow(tree$edge)){
    A<-pp$yy[tree$edge[i,1]]
    K<-pp$yy[tree$edge[i,2]]
    if(i==1) dy<-abs(A-K)
    B<-b*Ntip(tree)/h
    t<-seq(pp$xx[tree$edge[i,1]],pp$xx[tree$edge[i,2]],
           length.out=res)
    t<-sort(c(t,t[1]+cumsum(tree$maps[[i]])))
    dd<-diff(range(t))
    M<-t[1] + if(m1*h>(m2*dd)) m2*dd else m1*h
    C<-1
    Yt<-c(A,sigmoid(t),K)
    t<-c(t[1],t,t[length(t)])
    COLS<-vector()
    bb<-setNames(t[1]+cumsum(tree$maps[[i]]),names(tree$maps[[i]]))
    for(j in 1:length(t))
      COLS[j]<-colors[names(bb[which(t[j]<=bb)])[1]]
    nn<-length(t)
    segments(t[1:(nn-1)],Yt[1:(nn-1)],x1=t[2:nn],y1=Yt[2:nn],
             col=COLS,lwd=lwd)
  }
  nulo<-dev.flush()
}




library(nationalparkcolors)
pal <- park_palette("Saguaro")
library(phytools)
packageVersion("phytools")
Data <- read.csv("min20516.csv")
Data<-read.csv(file="mamPhyloCut.csv", row.names = 10)
Data<- read.csv("mamPhylocut.csv")
Data<- filter(Data, Clade == "Mammalia")
tree <- read.tree("min20Fixed516.nwk")

Data$Species <- gsub(" ", "_", Data$Species)
includedSpecies <- Data$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$Species %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


NeoplasiaPrevalence <- setNames(Data$NeoplasiaPrevalence, rownames(Data))
mamContMap <- contMap(pruned.tree, NeoplasiaPrevalence, plot = F)
mamContMap <- setMap(mamContMap, invert = T)
sigmoidPhylogram(pruned.tree,fsize=0.3,ftype="i",lwd=7,
                 ylim=c(-5,Ntip(pruned.tree)))

sigmoidPhylogram(mamContMap$tree,colors=mamContMap$cols,
                 ftype="off",lwd=5,
                 xlim=get("last_plot.phylo",envir=.PlotPhyloEnv)$x.lim,
                 ylim=c(-5,Ntip(pruned.tree)),
                 add=TRUE)
