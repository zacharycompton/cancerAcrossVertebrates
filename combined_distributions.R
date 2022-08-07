#open libraries
library(patchwork)
library(tidyverse)
library(cowplot)


#read csv
AgeRisk <- read.csv("cleanPath.min20.062822.csv")
#filter for sauropsids
AgeRisk_Sauropsida <- filter(AgeRisk, is.element(Clade, c("Sauropsida")))
AgeRisk_Sauropsida <- filter(AgeRisk_Sauropsida,proportion_lifespan >= 0)
#filter for mammals
AgeRisk_Mammals<- filter(AgeRisk, is.element(Class, c("Mammalia")))
AgeRisk_Mammals <- filter(AgeRisk_Mammals, proportion_lifespan>= 0)
#filter for amphibians
AgeRisk_Amphibia <- filter(AgeRisk, is.element(Class, c("Amphibia")))
AgeRisk_Amphibia <- filter(AgeRisk_Amphibia, proportion_lifespan >= 0)


#Density plot for percentage of age of death for Sauropsids

saur <- ggplot(AgeRisk_Sauropsida, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Sauropsida") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "none")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

#Density plot for Mammals
mam <- ggplot(AgeRisk_Mammals, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Mammalia") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "none")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

#Density plot for Amphibians
amph <- ggplot(AgeRisk_Amphibia, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Amphibia") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "none")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

#create a layout matrix
layout<-"
#AA#
BBCC"

#plot all plots onto one image
combo<- mam+saur+amph+plot_layout(design=layout,guides="collect", widths=c(10,10,10)) & theme(legend.position = "bottom")
combo+plot_annotation(caption="Lifespan Data obtained from panTHERIA,which often estimates wild populations\nY axis is probability density")

