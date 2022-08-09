library(tidyverse)
library(phytools)
library(OUwie)
library(ape)
library(ggpubr)
library(gridGraphics)
library(patchwork)
library(cowplot)

## read data and tree frosm file
##Mammals


Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
Data <- Data[,c(9,10,17,13),drop=FALSE] 
tree<-read.tree(file="commontree.nwk")
Data[Data < 0] <-NA
Data <- na.omit(Data)

Data$common_name <- gsub(" ", "_", Data$common_name)
includedSpecies <- Data$common_name
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
Data$Keep <- Data$common_name %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


rownames(Data)<-Data$common_name
Data <- Data[,c(3,4),drop=FALSE] 
name.check(pruned.tree,Data)


colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)

par(mfrow=c(3,1))

## extract scores for first two PCs
anole.pc<-(matData)[,1:2]
## set background black & foreground to white for plotting
par(fg="white",bg="black")
## compute max tree height and number of traits
h<-max(nodeHeights(pruned.tree))
m<-ncol(anole.pc)
## set colors
cols<-c("#631879FF","#008b45ff")
## set x & y limits
xlim<-ylim<-1.2*c(-h,h)+c(-1,1)*0.15*m*h+0.2*c(-h,h)
ylim<-c(0,ylim[2])
## plot tree
plotTree(pruned.tree,type="fan",xlim=xlim,ylim=ylim,
         lwd=1,ftype="i",fsize=0.5,part=0.5,color="white")
title(main = "A",col.main= "white",adj = .1, line = -1)

## add traits
for(i in 1:m){
  tt<-pruned.tree ## copy tree
  ## extend tip edges
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+
    0.35*h+(i-1)*0.2*h
  ## plot transparently
  plotTree(tt,color="transparent",type="fan",xlim=xlim,
           ylim=ylim,lwd=1,ftype="off",add=TRUE,part=0.5)
  pp1<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## extend again
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+
    0.2*h
  ## clip plot area
  clip(x1=par()$usr[1],x2=par()$usr[2],y1=ylim[1],y2=par()$usr[4])
  ## draw circle
  plotrix::draw.circle(0,0,radius=h+0.35*h+(i-1)*0.2*h,
                       border="#B8B8B8")
  ## unclip
  clip(x1=par()$usr[1],x2=par()$usr[2],y1=par()$usr[3],
       y2=par()$usr[4])
  ## plot transparently again
  plotTree(tt,color="transparent",type="fan",xlim=xlim,
           ylim=ylim,lwd=1,ftype="off",add=TRUE,part=0.5)
  pp2<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## graph lines
  par(lend=1)
  for(j in 1:Ntip(pruned.tree)){
    ii<-which(rownames(anole.pc)==pruned.tree$tip.label[j])
    dx<-(pp2$xx[j]-pp1$xx[j])*anole.pc[ii,i]
    dy<-(pp2$yy[j]-pp1$yy[j])*anole.pc[ii,i]
    lines(pp1$xx[j]+c(0,dx),pp1$yy[j]+c(0,dy),lwd=8,
          col=cols[i])
  }
}
## create custom legend
title(main="A", col = "white")
xx<-rep(0.85*par()$usr[4],m)
yy<-0.95*par()$usr[4]-1.5*1:m*strheight("W")
text(xx,yy,colnames(anole.pc),pos=4,cex=0.8)
#scale.bar<-sqrt(diag(anole.pca$Eval)[1:3])*0.2*h
#xx2<-xx-scale.bar
segments(x0=xx,y0=yy,x1=200,y1=yy,lwd=8,col=cols)
text(x=xx,y=0.95*par()$usr[4],"Factors",
     pos=4,cex=0.8)



## read data and tree frosm file
##Sauropsids


Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
Data <- Data[,c(9,10,17,13),drop=FALSE] 
tree<-read.tree(file="commontree.nwk")
Data[Data < 0] <-NA
Data <- na.omit(Data)

Data$common_name <- gsub(" ", "_", Data$common_name)
includedSpecies <- Data$common_name
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
Data$Keep <- Data$common_name %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


rownames(Data)<-Data$common_name
Data <- Data[,c(3,4),drop=FALSE] 
name.check(pruned.tree,Data)



colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)

anole.pc<-(matData)[,1:2]
## set background black & foreground to white for plotting
par(fg="white",bg="black")
## compute max tree height and number of traits
h<-max(nodeHeights(pruned.tree))
m<-ncol(anole.pc)
## set colors
 cols<-c("#631879FF","#008b45ff")
## set x & y limits
xlim<-ylim<-1.2*c(-h,h)+c(-1,1)*0.15*m*h+0.2*c(-h,h)
ylim<-c(0,ylim[2])
## plot tree
plotTree(pruned.tree,type="fan",xlim=xlim,ylim=ylim,
         lwd=1,ftype="i",fsize=0.5,part=0.5,color="white")
title(main = "B",col.main= "white",adj = .1, line = -1)
## add traits
for(i in 1:m){
  tt<-pruned.tree ## copy tree
  ## extend tip edges
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+
    0.35*h+(i-1)*0.2*h
  ## plot transparently
  plotTree(tt,color="transparent",type="fan",xlim=xlim,
           ylim=ylim,lwd=1,ftype="off",add=TRUE,part=0.5)
  pp1<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## extend again
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+
    0.2*h
  ## clip plot area
  clip(x1=par()$usr[1],x2=par()$usr[2],y1=ylim[1],y2=par()$usr[4])
  ## draw circle
  plotrix::draw.circle(0,0,radius=h+0.35*h+(i-1)*0.2*h,
                       border="#B8B8B8")
  ## unclip
  clip(x1=par()$usr[1],x2=par()$usr[2],y1=par()$usr[3],
       y2=par()$usr[4])
  ## plot transparently again
  plotTree(tt,color="transparent",type="fan",xlim=xlim,
           ylim=ylim,lwd=1,ftype="off",add=TRUE,part=0.5)
  pp2<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## graph lines
  par(lend=1)
  for(j in 1:Ntip(pruned.tree)){
    ii<-which(rownames(anole.pc)==pruned.tree$tip.label[j])
    dx<-(pp2$xx[j]-pp1$xx[j])*anole.pc[ii,i]
    dy<-(pp2$yy[j]-pp1$yy[j])*anole.pc[ii,i]
    lines(pp1$xx[j]+c(0,dx),pp1$yy[j]+c(0,dy),lwd=8,
          col=cols[i])
  }
}
## legend for individual tree
#xx<-rep(0.85*par()$usr[4],m)
#yy<-0.95*par()$usr[4]-1.5*1:m*strheight("W")
#text(xx,yy,colnames(anole.pc),pos=4,cex=0.8)
#segments(x0=xx,y0=yy,x1=350,y1=yy,lwd=8,col=cols)
#text(x=xx,y=0.95*par()$usr[4],"Factors",
#     pos=4,cex=0.8)


## read data and tree frosm file
##Amphibians

Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Amphibia")))
Data <- Data[,c(9,10,17,13),drop=FALSE] 
tree<-read.tree(file="amphtree.nwk")
Data[Data < 0] <-NA
Data <- na.omit(Data)

Data$common_name <- gsub(" ", "_", Data$common_name)
Data$common_name <- gsub("'", "", Data$common_name)
rownames(Data)<-Data$common_name
Data <- Data[,c(3,4),drop=FALSE] 
name.check(tree,Data)


colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)


## run phylogenetic PCA
#anole.pca<-phyl.pca(pruned.tree,matData,mode="corr")
#anole.pca
## extract scores for first two PCs
anole.pc<-(matData)[,1:2]
## normalize to have the same variance
#anole.norm<-anole.pc/matrix(rep(apply(anole.pc,2,sd),nrow(anole.pc)),
#                            nrow(anole.pc),ncol(anole.pc),byrow=TRUE)*0.2
## set background black & foreground to white for plotting
par(fg="white",bg="black")
## compute max tree height and number of traits
h<-max(nodeHeights(tree))
m<-ncol(anole.pc)
## set colors
 cols<-c("#631879FF","#008b45ff")
## set x & y limits
xlim<-ylim<-1.2*c(-h,h)+c(-1,1)*0.15*m*h+0.2*c(-h,h)
ylim<-c(0,ylim[2])
## plot tree
plotTree(tree,type="fan",xlim=xlim,ylim=ylim,
         lwd=1,ftype="i",fsize=0.5,part=0.5,color="white")
title(main = "C",col.main= "white",adj = .1, line = -1)
## add traits
for(i in 1:m){
  tt<-tree ## copy tree
  ## extend tip edges
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+
    0.35*h+(i-1)*0.2*h
  ## plot transparently
  plotTree(tt,color="transparent",type="fan",xlim=xlim,
           ylim=ylim,lwd=1,ftype="off",add=TRUE,part=0.5)
  pp1<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## extend again
  tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]<-
    tt$edge.length[which(tt$edge[,2]<=Ntip(tt))]+
    0.2*h
  ## clip plot area
  clip(x1=par()$usr[1],x2=par()$usr[2],y1=ylim[1],y2=par()$usr[4])
  ## draw circle
  plotrix::draw.circle(0,0,radius=h+0.35*h+(i-1)*0.2*h,
                       border="#B8B8B8")
  ## unclip
  clip(x1=par()$usr[1],x2=par()$usr[2],y1=par()$usr[3],
       y2=par()$usr[4])
  ## plot transparently again
  plotTree(tt,color="transparent",type="fan",xlim=xlim,
           ylim=ylim,lwd=1,ftype="off",add=TRUE,part=0.5)
  pp2<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  ## graph lines
  par(lend=1)
  for(j in 1:Ntip(tree)){
    ii<-which(rownames(anole.pc)==tree$tip.label[j])
    dx<-(pp2$xx[j]-pp1$xx[j])*anole.pc[ii,i]
    dy<-(pp2$yy[j]-pp1$yy[j])*anole.pc[ii,i]
    lines(pp1$xx[j]+c(0,dx),pp1$yy[j]+c(0,dy),lwd=8,
          col=cols[i])
  }
}
## legend for individual tree
#xx<-rep(0.85*par()$usr[4],m)
#yy<-0.95*par()$usr[4]-1.5*1:m*strheight("W")
#text(xx,yy,colnames(anole.pc),pos=4,cex=0.8)
#segments(x0=xx,y0=yy,x1=350,y1=yy,lwd=8,col=cols)
#text(x=xx,y=0.95*par()$usr[4],"Factors",
  #   pos=4,cex=0.8)



