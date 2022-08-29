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


Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Mammalia")))
Data <- Data[,c(6,9,10,17,13),drop=FALSE] 
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

specData<-Data

rownames(Data)<-Data$common_name
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)

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


colors<-c("#1B1919FF","#3B4992FF","#EE0000FF","#008280FF","#A20056FF","#5F559BFF")

names(colors)<-0:5

pruned.tree<-paintSubTree(pruned.tree,node=110,state="1",anc="0")
pruned.tree<-paintSubTree(pruned.tree,node=91,state="2")
pruned.tree<-paintSubTree(pruned.tree,node=78,state="3")
pruned.tree<-paintSubTree(pruned.tree,node=118,state="4")
pruned.tree<-paintSubTree(pruned.tree,node=68,state="5")


plotSimmap(pruned.tree,colors,type="fan",xlim=xlim,ylim=ylim,
           lwd=1,ftype="i",fsize=1,part=0.5)

par(fg="#3B4992FF")
arc.cladelabels(text="Rodentia",cex = .8,node=findMRCA(pruned.tree, rodentia$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

par(fg="#EE0000FF")
arc.cladelabels(text="Artiodactyla",cex = .8,node=findMRCA(pruned.tree, artio$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Carnivora",cex = .8,node=findMRCA(pruned.tree, carnivora$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

#arc.cladelabels(text="Cetacea",cex = .8,node=findMRCA(pruned.tree, cetacea$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text=" Chiroptera",cex = .8,node=findMRCA(pruned.tree, chiroptera$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Didelphimorphia",cex=.5,node=findMRCA(pruned.tree, didelphimorphia$common_name),ln.offset=1.7,lab.offset=1.8, orientation = "horizontal")
#arc.cladelabels(text="Eulipotyphla",node=which(pruned.tree$tip.label=="Four-toed_hedgehog"), orientation="horizontal",ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Hyracoidea",fsize=.4,node=65,orientation="horizontal",ln.offset=1.45,lab.offset=1.5,mark.node=FALSE)
#arc.cladelabels(text="Lagomorpha",cex.sub=.1,node=which(pruned.tree$tip.label=="Domestic_rabbit"),orientation="horizontal",ln.offset=1.45,lab.offset=1.45, mark.node=FALSE)
#arc.cladelabels(text="Perissodactyla",cex = .8,node=which(pruned.tree$tip.label=="Grevys_zebra"),orientation="horizontal",ln.offset=1.45,lab.offset=1.5, mark.node=FALSE)
par(fg="#A20056FF")
arc.cladelabels(text="Primates",cex = .8,node=findMRCA(pruned.tree, primates$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

#arc.cladelabels(text="Proboscidea",cex = .8,node=which(pruned.tree$tip.label=="Asian_elephant"),orientation="horizontal",ln.offset=1.45,lab.offset=1.5, mark.node=FALSE)
par(fg="#5F559BFF")
arc.cladelabels(text="Diprotodontia",cex = .8,node=findMRCA(pruned.tree, diprotodontia$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

#arc.cladelabels(text="Cingulata",node=which(pruned.tree$tip.label=="Nine-banded_armadillo"),ln.offset=1.45,lab.offset=1.5,fsize=0.5, orientation = "horizontal",mark.node=FALSE)

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


Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Sauropsida")))
Data <- Data[,c(6,9,10,17,13),drop=FALSE] 
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

specData<-Data




rownames(Data)<-Data$common_name
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)



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


plotTree(pruned.tree,type="fan",xlim=xlim,ylim=ylim,
         lwd=1,ftype="i",fsize=1,part=0.5)
colors<-c("#1B1919FF","#3B4992FF","#EE0000FF","#BB0021FF","#008280FF","#A20056FF","#5F559BFF")

pruned.tree<-paintSubTree(pruned.tree,node=165,state="1",anc="0")
pruned.tree<-paintSubTree(pruned.tree,node=174,state="2")
pruned.tree<-paintSubTree(pruned.tree,node=157,state="3")
pruned.tree<-paintSubTree(pruned.tree,node=136,state="4")
pruned.tree<-paintSubTree(pruned.tree,node=130,state="5")
pruned.tree<-paintSubTree(pruned.tree,node=92,state="6")

names(colors)<-0:6

plotSimmap(pruned.tree,colors,type="fan",xlim=xlim,ylim=ylim,
           lwd=1,ftype="i",fsize=1,part=0.5)


#arc.cladelabels(text="Accipitriformes",cex = .8,node=findMRCA(pruned.tree, accipitriformes$common_name),ln.offset=1.7,lab.offset=1.8)
par(fg="#3B4992FF")
arc.cladelabels(text="Anseriformes",cex = .8,node=findMRCA(pruned.tree, anseriformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)

#arc.cladelabels(text="Bucerotiformes",cex = .8,node=findMRCA(pruned.tree, bucerotiformesiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Charadriiformes",cex = .8,node=findMRCA(pruned.tree, charadriiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Caprimulgiformes",cex = .8,node=findMRCA(pruned.tree, caprimulgiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text=" Ciconiiformes",cex = .8,node=findMRCA(pruned.tree, ciconiiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="    Columbiformes",cex = .8,node=findMRCA(pruned.tree, columbiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Crocodilia",cex = .8,node=findMRCA(pruned.tree, crocodilia$common_name),ln.offset=1.7,lab.offset=1.75)
par(fg="#EE0000FF")
arc.cladelabels(text="Galliformes",cex = .8,node=findMRCA(pruned.tree, galliformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)
par(fg="#BB0021FF")
arc.cladelabels(text="Pelecaniformes",cex = .8,node=findMRCA(pruned.tree, pelecaniformes$common_name),ln.offset=1.7,lab.offset=1.75, mark.node=FALSE,mark.node=FALSE)

par(fg="#008280FF")
arc.cladelabels(text="Passeriformes",cex = .8,node=findMRCA(pruned.tree, passeriformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)
#arc.cladelabels(text="piciformes",cex = .8,node=findMRCA(pruned.tree, piciformes$common_name),ln.offset=1.7,lab.offset=1.75)
par(fg="#A20056FF")
arc.cladelabels(text="Psittaciformes",cex = .8,node=findMRCA(pruned.tree, psittaciformes$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE)
#arc.cladelabels(text="Rheiformes",cex = .8,node=findMRCA(pruned.tree, rheiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Sphenisciformes",cex = .8,node=findMRCA(pruned.tree, sphenisciformes$common_name),ln.offset=1.7,lab.offset=1.75)
par(fg="#5F559BFF")
arc.cladelabels(text="Squamata",cex = .8,node=findMRCA(pruned.tree, squamata$common_name),ln.offset=1.7,lab.offset=1.75,mark.node=FALSE,mark.node=FALSE)
#arc.cladelabels(text="Strigiformes",cex = .8,node=findMRCA(pruned.tree, strigiformes$common_name),ln.offset=1.7,lab.offset=1.75)
#arc.cladelabels(text="Testudines",cex = .8,node=findMRCA(pruned.tree, testudines$common_name),ln.offset=1.7,lab.offset=1.75)

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

Data<-read.csv(file="min20516.csv")
Data<- filter(Data, is.element(Clade, c("Amphibia")))
Data <- Data[,c(6,9,10,17,13),drop=FALSE] 
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

specData<-Data

rownames(Data)<-Data$common_name
Data <- Data[,c(4,5),drop=FALSE] 
name.check(pruned.tree,Data)


colnames(Data) <- c('Malignancy Prevalence', 'Neoplasia Prevalence')

matData<-as.matrix(Data)


## run phylogenetic PCA
#anole.pca<-phyl.pca(pruned.pruned.tree,matData,mode="corr")
#anole.pca
## extract scores for first two PCs
anole.pc<-(matData)[,1:2]
## normalize to have the same variance
#anole.norm<-anole.pc/matrix(rep(apply(anole.pc,2,sd),nrow(anole.pc)),
#                            nrow(anole.pc),ncol(anole.pc),byrow=TRUE)*0.2
## set background black & foreground to white for plotting
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


colors<-c("#1B1919FF","#3B4992FF","#EE0000FF")

pruned.tree<-paintSubTree(pruned.tree,node=20,state="1",anc="0")
pruned.tree<-paintSubTree(pruned.tree,node=17,state="2")

names(colors)<-0:2


plotSimmap(pruned.tree,colors,type="fan",xlim=xlim,ylim=ylim,
           lwd=1,ftype="i",fsize=1,part=0.5)

par(fg="#3B4992FF")
arc.cladelabels(text="Anura",cex = .8,node=findMRCA(pruned.tree, anura$common_name),ln.offset=1.7,lab.offset=1.75, mark.node=FALSE)
par(fg="#EE0000FF")
arc.cladelabels(text="Caudata",cex = .8,node=findMRCA(pruned.tree, caudata$common_name),ln.offset=1.7,lab.offset=1.75, mark.node=FALSE)
#arc.cladelabels(text="Gymnophiona",cex.sub=.1,node=which(pruned.tree$tip.label=="Gaboon_caecilian"),orientation="horizontal",ln.offset=1.45,lab.offset=1.45, mark.node=FALSE)

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



