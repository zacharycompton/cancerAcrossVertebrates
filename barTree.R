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


Data<-read.csv(file="min20-2022.05.16.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
Data <- Data[,c(6,9,10,17,13),drop=FALSE] 
Data[Data < 0] <-NA
Data <- na.omit(Data)

tree<-read.tree(file="commontree.nwk")

Data$common_name <- gsub(" ", "_", Data$common_name)
includedSpecies <- Data$common_name
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$common_name %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]

specData<-Data

rownames(Data)<-Data$common_name
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)


#create malignancy and neoplasia matrix

colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)

#par(mfrow=c(3,1))

## extract scores for first two PCs
anole.pc<-(matData)[,1:2]
## set background black & foreground to white for plotting
par(fg="black",bg="white")
## compute max tree height and number of traits
h<-max(nodeHeights(pruned.tree))
m<-ncol(anole.pc)
## set colors
cols<-c("#631879FF","#008b45ff")
## set x & y limits
xlim<-ylim<-1.2*c(-h,h)+c(-1,1)*0.15*m*h+0.2*c(-h,h)
ylim<-c(0,ylim[2])
## plot tree
#subset for labeling
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
#label function
arc.cladelabels<-function(tree=NULL,text,node,ln.offset=1.02,
                          lab.offset=1.06,cex=1,orientation="curved",...){
  obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  if(obj$type!="fan") stop("method works only for type=\"fan\"")
  h<-max(sqrt(obj$xx^2+obj$yy^2))
  if(hasArg(mark.node)) mark.node<-list(...)$mark.node
  else mark.node<-TRUE
  if(mark.node) points(obj$xx[node],obj$yy[node],pch=21,
                       bg="red")
  if(is.null(tree)){
    tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
               Nnode=obj$Nnode)
    class(tree)<-"phylo"
  }
  d<-getDescendants(tree,node)
  d<-sort(d[d<=Ntip(tree)])
  deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
  ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
  deg[ii]<-180+deg[ii]
  ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
  deg[ii]<-360+deg[ii]
  draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),
           deg2=max(deg))
  if(orientation=="curved")
    arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex)
  else if(orientation=="horizontal"){
    x0<-lab.offset*cos(median(deg)*pi/180)*h
    y0<-lab.offset*sin(median(deg)*pi/180)*h
    text(x=x0,y=y0,label=text,
         adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
         offset=0)
  }
}
#plot with labels


#create color list

colors<-c("#1B1919FF","#3B4992FF","#EE0000FF","#008280FF","#BB0021FF","#A20056FF","#5F559BFF")

names(colors)<-0:6

#color code tree branches

pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, rodentia$common_name),state="1",anc="0")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, artio$common_name),state="2")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, carnivora$common_name),state="3")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, chiroptera$common_name),state="4")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, primates$common_name),state="5")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, diprotodontia$common_name),state="6")


plotSimmap(pruned.tree,colors,type="fan",xlim=xlim,ylim=ylim,
           lwd=1,ftype="i",fsize=1,part=0.5)


#order labels

par(fg="#3B4992FF")
arc.cladelabels(text="Rodentia",cex = .8,node=findMRCA(pruned.tree, rodentia$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

par(fg="#EE0000FF")
arc.cladelabels(text="Artiodactyla",cex = .8,node=findMRCA(pruned.tree, artio$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Carnivora",cex = .8,node=findMRCA(pruned.tree, carnivora$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)


par(fg="#BB0021FF")
arc.cladelabels(text=" Chiroptera",cex = .8,node=findMRCA(pruned.tree, chiroptera$common_name),ln.offset=1.7,lab.offset=1.75,mark.node =FALSE)

par(fg="#A20056FF")
arc.cladelabels(text="Primates",cex = .8,node=findMRCA(pruned.tree, primates$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

par(fg="#5F559BFF")
arc.cladelabels(text="Diprotodontia",cex = .8,node=findMRCA(pruned.tree, diprotodontia$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)


title(main = "A",col.main= "black",adj = .1, line = -1)


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
xx<-rep(0.85*par()$usr[4],m)
yy<-0.95*par()$usr[4]-1.5*1:m*strheight("W")
text(xx+20,yy,colnames(anole.pc),pos=4,cex=0.8)
#scale.bar<-sqrt(diag(anole.pca$Eval)[1:3])*0.2*h
#xx2<-xx-scale.bar
segments(x0=xx,y0=yy,x1=260,y1=yy,lwd=8,col=cols)





## read data and tree frosm file
##Sauropsids


Data<-read.csv(file="min20-2022.05.16.csv")
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
Data <- Data[,c(4,6,9,10,17,13),drop=FALSE] 
Data[Data < 0] <-NA
Data <- na.omit(Data)

tree<-read.tree(file="commontree.nwk")

Data$common_name <- gsub(" ", "_", Data$common_name)
includedSpecies <- Data$common_name
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$common_name %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]

specData<-Data

rownames(Data)<-Data$common_name
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)


#create malignancy and neoplasia matrix

colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)

anole.pc<-(matData)[,1:2]
## set background black & foreground to white for plotting
par(fg="black",bg="white")
## compute max tree height and number of traits
h<-max(nodeHeights(pruned.tree))
m<-ncol(anole.pc)
## set colors
cols<-c("#631879FF","#008b45ff")
## set x & y limits
xlim<-ylim<-1.2*c(-h,h)+c(-1,1)*0.15*m*h+0.2*c(-h,h)
ylim<-c(0,ylim[2])
## plot tree


#filter by order

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
aves<-filter(specData, is.element(Class, c("Aves")))


plotTree(pruned.tree,type="fan",xlim=xlim,ylim=ylim,
         lwd=1,ftype="i",fsize=1,part=0.5)


#color tree branches

colors<-c("#1B1919FF","#008280FF","#EE0000FF","#BB0021FF","#3B4992FF","#A20056FF","#5F559BFF")

pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, aves$common_name),state="1",anc="0")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, galliformes$common_name),state="2")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, pelecaniformes$common_name),state="3")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, passeriformes$common_name),state="4")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, psittaciformes$common_name),state="5")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, squamata$common_name),state="6")

names(colors)<-0:6

plotSimmap(pruned.tree,colors,type="fan",xlim=xlim,ylim=ylim,
           lwd=1,ftype="i",fsize=1,part=0.5)


#add labels

par(fg="#008280FF")
arc.cladelabels(text="Aves",cex = .8,node=findMRCA(pruned.tree, aves$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

# par(fg="#EE0000FF")
# arc.cladelabels(text="Galliformes",cex = .8,node=findMRCA(pruned.tree, galliformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)
# 
# par(fg="#BB0021FF")
# arc.cladelabels(text="Pelecaniformes",cex = .8,node=findMRCA(pruned.tree, pelecaniformes$common_name),ln.offset=1.7,lab.offset=1.75, mark.node=FALSE,mark.node=FALSE)
# 
# par(fg="#008280FF")
# arc.cladelabels(text="Passeriformes",cex = .8,node=findMRCA(pruned.tree, passeriformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)
# 
# par(fg="#A20056FF")
# arc.cladelabels(text="Psittaciformes",cex = .8,node=findMRCA(pruned.tree, psittaciformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

par(fg="#5F559BFF")
arc.cladelabels(text="Squamata",cex = .8,node=findMRCA(pruned.tree, squamata$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE,mark.node=FALSE)

title(main = "B",col.main= "black",adj = 0, line = -1)


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

Data<-read.csv(file="min20-2022.05.16.csv")
Data<- filter(Data, is.element(Clade, c("Amphibia")))
Data <- Data[,c(6,9,10,17,13),drop=FALSE] 
Data[Data < 0] <-NA
Data <- na.omit(Data)

tree<-read.tree(file="commontree.nwk")

Data$common_name <- gsub(" ", "_", Data$common_name)
includedSpecies <- Data$common_name
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$common_name %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]

specData<-Data

rownames(Data)<-Data$common_name
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)

#create malignancy and neoplasia matrix
colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)


anole.pc<-(matData)[,1:2]

par(fg="black",bg="white")
## compute max pruned.tree height and number of traits
h<-max(nodeHeights(pruned.tree))
m<-ncol(anole.pc)
## set colors
cols<-c("#631879FF","#008b45ff")
## set x & y limits
xlim<-ylim<-1.2*c(-h,h)+c(-1,1)*0.15*m*h+0.2*c(-h,h)
ylim<-c(0,ylim[2])



## plot pruned.tree
plotTree(pruned.tree,type="fan",xlim=xlim,ylim=ylim,
         lwd=1,ftype="i",fsize=1,part=0.5,color="black")

anura<-filter(specData, is.element(Orders, c("Anura")))
caudata<-filter(specData, is.element(Orders, c("Caudata")))
gymnophiona<-filter(specData, is.element(Orders, c("Gymnophiona")))

#create color list
colors<-c("#1B1919FF","#3B4992FF","#EE0000FF")


#color tree branches
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, anura$common_name),state="1",anc="0")
pruned.tree<-paintSubTree(pruned.tree,node=findMRCA(pruned.tree, caudata$common_name),state="2")

names(colors)<-0:2


plotSimmap(pruned.tree,colors,type="fan",xlim=xlim,ylim=ylim,
           lwd=1,ftype="i",fsize=1,part=0.5)

par(fg="#3B4992FF")
arc.cladelabels(text="Anura",cex = .8,node=findMRCA(pruned.tree, anura$common_name),ln.offset=1.7,lab.offset=1.75, mark.node=FALSE)
par(fg="#EE0000FF")
arc.cladelabels(text="Caudata",cex = .8,node=findMRCA(pruned.tree, caudata$common_name),ln.offset=1.7,lab.offset=1.75, mark.node=FALSE)

title(main = "C",col.main= "black",adj = .1, line = -1)
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
#   pos=4,cex=0.8)



