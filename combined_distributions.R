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
  ggtitle("Sauropsids Neoplasia") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))
ggsave("saurage.png",width=13, height=10, limitsize=FALSE,bg="white")

#Density plot for Mammals
mam <- ggplot(AgeRisk_Mammals, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Mammals Neoplasia ") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))
ggsave("mammage.png",width=13, height=10, limitsize=FALSE,bg="white")


#Density plot for Amphibians
amph <- ggplot(AgeRisk_Amphibia, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Amphibians Neoplasia") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

ggsave("amphage.png",width=13, height=10, limitsize=FALSE,bg="white")

#create a layout matrix
layout<-"
AA
BB
CC"

#plot all plots onto one image
combo<- mam+saur+amph+plot_layout(design=layout,guides="collect", widths=c(10,10,10)) & theme(legend.position = "bottom")
combo+plot_annotation(caption="Lifespan Data obtained from panTHERIA,which often estimates wild populations\nY axis is probability density")

ggsave("density.pdf",width=5.5, height=9, limitsize=FALSE,bg="white")




#Density plot for percentage of age of death for Sauropsids
AgeRisk <- AgeRisk[,c(14,19,51,52),drop=FALSE] 

#used for sanity check
#newAgeRisk<-AgeRisk %>% rowwise() %>%
 # mutate(malandneo = sum(c_across(Masspresent:Malignant)))

newAgeRisk<-mutate(AgeRisk, malandneo = ifelse(Malignant == 1, 1,
                                  ifelse(Malignant == -1|Malignant == 0, 0,NA)))

view(newAgeRisk)



AgeRisk_Sauropsida <- filter(newAgeRisk, is.element(Clade, c("Sauropsida")))
AgeRisk_Sauropsida <- filter(AgeRisk_Sauropsida,proportion_lifespan >= 0)
#filter for mammals
AgeRisk_Mammals<- filter(newAgeRisk, is.element(Clade, c("Mammalia")))
AgeRisk_Mammals <- filter(AgeRisk_Mammals, proportion_lifespan>= 0)
#filter for amphibians
AgeRisk_Amphibia <- filter(newAgeRisk, is.element(Clade, c("Amphibia")))
AgeRisk_Amphibia <- filter(AgeRisk_Amphibia, proportion_lifespan >= 0)

saurMal <- ggplot(AgeRisk_Sauropsida, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((malandneo)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Sauropsids Malignancy") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))
#ggsave("saurage.png",width=13, height=10, limitsize=FALSE,bg="white")

#Density plot for Mammals
mamMal <- ggplot(AgeRisk_Mammals, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((malandneo)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Mammals Malignancy") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))
#ggsave("mammage.png",width=13, height=10, limitsize=FALSE,bg="white")


#Density plot for Amphibians
amphMal <- ggplot(AgeRisk_Amphibia, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((malandneo)))) +
  geom_density(alpha=0.25) + 
  ggtitle("Amphibians Malignancy") + 
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

#ggsave("amphage.png",width=13, height=10, limitsize=FALSE,bg="white")

#create a layout matrix
layout<-"
AADD
BBEE
CCFF"

#plot all plots onto one image
combo<- mam+saur+amph+mamMal+saurMal+amphMal+plot_layout(design=layout,guides="collect", widths=c(10,10,10)) & theme(legend.position = "bottom")
combo+plot_annotation(caption="Lifespan Data obtained from panTHERIA,which often estimates wild populations\nY axis is probability density")

ggsave("Maldensity.pdf",width=11, height=9, limitsize=FALSE,bg="white")



 

