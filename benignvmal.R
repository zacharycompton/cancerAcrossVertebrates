#open libraries
library(patchwork)
library(tidyverse)
library(cowplot)


#read csv
cleanPath <- read.csv("cleanPath.min20.062822.csv")
min20Data <- read.csv("min20516.csv")
view(min20Data)
malben<-lm(MalignancyPrevalence~BenignPrevalence,data = min20Data) 
summary(malben)
summary(malben)$r.squared
summary(malben)$coefficients 

benign<-filter(cleanPath, is.element(Malignant, c(0)))
malignant<-filter(cleanPath, is.element(Malignant, c(1)))
ks.test(benign$proportion_lifespan,malignanct$proportion_lifespan)
lm(benign$Neop )
cleanPath<- filter(cleanPath,Malignant >= 0)






bvm<- ggplot(cleanPath, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Malignant)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Slow Life History",
          subtitle="All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c( "No Tumor Found","Benign", "Malignant"))

bvm