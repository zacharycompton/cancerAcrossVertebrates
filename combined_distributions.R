install.packages('gapminder')
install.packages("patchwork")
library(gapminder)
library(patchwork)
library(tidyverse)
library(gridExtra)
library(cowplot)
library(ggpubr)

AgeRisk <- read.csv("cleanPath.min20.062822.csv")
AgeRisk_Sauropsida <- filter(AgeRisk, is.element(Clade, c("Sauropsida")))
AgeRisk_Sauropsida <- filter(AgeRisk_Sauropsida,proportion_lifespan >= 0)
AgeRisk_Mammals<- filter(AgeRisk, is.element(Class, c("Mammalia")))
AgeRisk_Mammals <- filter(AgeRisk_Mammals, proportion_lifespan>= 0)
AgeRisk_Amphibia <- filter(AgeRisk, is.element(Class, c("Amphibia")))
AgeRisk_Amphibia <- filter(AgeRisk_Amphibia, proportion_lifespan >= 0)


#Density plot of age(months) vs Mass Present for Reptilia

saur <- ggplot(AgeRisk_Sauropsida, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Sauropsida") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "none")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

mam <- ggplot(AgeRisk_Mammals, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Mammalia") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "none")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

amph <- ggplot(AgeRisk_Amphibia, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Amphibia") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "none")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))


layout<-"
#AA#
BBCC"

combo<- mam+saur+amph+plot_layout(design=layout,guides="collect", widths=c(10,10,10)) & theme(legend.position = "bottom")
combo+plot_annotation(caption="Lifespan Data obtained from panTHERIA,which often estimates wild populations\nY axis is probability density")

ggplot(AgeRisk_Sauropsida, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density()
