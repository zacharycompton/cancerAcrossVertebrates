#open libraries
library(patchwork)
library(tidyverse)
library(cowplot)


#read csv
fastSlow <- read.csv("fastslow.csv")
#filter for wgtlong
fastSlow <- filter(fastSlow,adult_weight >= 0)
fastSlow <- filter(fastSlow,max_longevity >= 0)

cutoff<-median(fastSlow$wgtxlong)


#fast == 1
#slow == 0
fastSlow$fastslow<-ifelse(fastSlow$wgtxlong>=cutoff,'0','1')

fast<-filter(fastSlow, is.element(fastslow, c(1)))
fast<-fast[fast$proportion_lifespan > 1, ]
slow<-filter(fastSlow, is.element(fastslow, c(0)))
slow<-slow[slow$proportion_lifespan > 1, ]
ks.test(fast$proportion_lifespan,slow$proportion_lifespan)

view(fastSlow)

fastvslowall<- ggplot(fastSlow, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density",
          subtitle = "All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))



fastvslowall

ggsave("fastvslowall.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#filter for mammals
fastSlowMamm <- filter(fastSlow,Clade == "Mammalia")

fastMamm<-filter(fastSlowMamm, is.element(fastslow, c(1)))
fastMamm<-fastMamm[fastMamm$proportion_lifespan > 1, ]
slowMamm<-filter(fastSlowMamm, is.element(fastslow, c(0)))
slowMamm<-slowMamm[slowMamm$proportion_lifespan > 1, ]
ks.test(fastMamm$proportion_lifespan,slowMamm$proportion_lifespan)

fastvslowMamm<- ggplot(fastSlowMamm, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density",
          subtitle = "Mammals") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

fastvslowMamm

ggsave("fastvslowMamm.pdf",width=13, height=10, limitsize=FALSE,bg="white")
#filer for sauropsids
fastSlowSaur <- filter(fastSlow,Clade == "Sauropsida")


fastSaur<-filter(fastSlowSaur, is.element(fastslow, c(1)))
fastSaur<-fastSaur[fastSaur$proportion_lifespan > 1, ]
slowSaur<-filter(fastSlowSaur, is.element(fastslow, c(0)))
slowSaur<-slowSaur[slowSaur$proportion_lifespan > 1, ]
ks.test(fastSaur$proportion_lifespan,slowSaur$proportion_lifespan)

fastvslowSaur<- ggplot(fastSlowSaur, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density",
          subtitle = "Sauropsids") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

fastvslowSaur

ggsave("fastvslowSaur.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#create a layout matrix
layout<-"
AA
BB
CC"

#plot all plots onto one image
combo<- fastvslowall+fastvslowMamm+fastvslowSaur+plot_layout(design=layout,guides="collect", widths=c(10,10,10)) & theme(legend.position = "bottom")
combo+plot_annotation(caption="Lifespan Data obtained from panTHERIA,which often estimates wild populations\nY axis is probability density")

ggsave("fastvslow.png",width=5.5, height=9, limitsize=FALSE,bg="white")

#fast tumor vs no tumor

fast <- filter(fastSlow,fastslow == 1)


notumorFoundFast<-filter(fast, is.element(Masspresent, c(0)))
notumorFoundFast<-notumorFoundFast[notumorFoundFast$proportion_lifespan > 1, ] 
tumorFoundFast<-filter(fast, is.element(Masspresent, c(1)))
tumorFoundFast<-tumorFoundFast[tumorFoundFast$proportion_lifespan > 1, ]
ks.test(notumorFoundFast$proportion_lifespan,tumorFoundFast$proportion_lifespan)

fasttumor <- ggplot(fast, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Fast Life History", 
          subtitle="All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

fasttumor

ggsave("fasttumor.pdf",width=13, height=10, limitsize=FALSE,bg="white")

#fast tumor vs no tumor

slow <- filter(fastSlow,fastslow == 0)

notumorFoundSlow<-filter(slow, is.element(Masspresent, c(0)))
notumorFoundSlow<-notumorFoundSlow[notumorFoundSlow$proportion_lifespan > 1, ]
tumorFoundSlow<-filter(slow, is.element(Masspresent, c(1)))
tumorFoundSlow<-tumorFoundSlow[tumorFoundSlow$proportion_lifespan > 1, ]
ks.test(notumorFoundSlow$proportion_lifespan,tumorFoundSlow$proportion_lifespan)


slowtumor<- ggplot(slow, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Slow Life History",
          subtitle="All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

slowtumor

ggsave("slowtumor.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#fast tumor vs no tumor mammals

fastMamm <- filter(fastSlowMamm,fastslow == 1)

notumorFoundFastMamm<-filter(fastMamm, is.element(Masspresent, c(0)))
notumorFoundFastMamm<-notumorFoundFastMamm[notumorFoundFastMamm$proportion_lifespan > 1, ]
tumorFoundFastMamm<-filter(fastMamm, is.element(Masspresent, c(1)))
tumorFoundFastMamm<-tumorFoundFastMamm[tumorFoundFastMamm$proportion_lifespan > 1, ]
ks.test(notumorFoundFastMamm$proportion_lifespan,tumorFoundFastMamm$proportion_lifespan)



fasttumorMamm <- ggplot(fastMamm, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Fast Life History", 
          subtitle="Mammals") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

fasttumorMamm

ggsave("fasttumorMamm.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#fast tumor vs no tumor

slowMamm <- filter(fastSlowMamm,fastslow == 0)

notumorFoundSlowMamm<-filter(slowMamm, is.element(Masspresent, c(0)))
notumorFoundSlowMamm<-notumorFoundSlowMamm[notumorFoundSlowMamm$proportion_lifespan > 1, ]
tumorFoundSlowMamm<-filter(slowMamm, is.element(Masspresent, c(1)))
tumorFoundSlowMamm<-tumorFoundSlowMamm[tumorFoundSlowMamm$proportion_lifespan > 1, ]
ks.test(notumorFoundSlowMamm$proportion_lifespan,tumorFoundSlowMamm$proportion_lifespan)


slowtumorMamm<- ggplot(slowMamm, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Slow Life History",
          subtitle="Mammals") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

slowtumorMamm

ggsave("slowtumorMamm.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#fast tumor vs no tumor Saurals

fastSaur <- filter(fastSlowSaur,fastslow == 1)

notumorFoundFastSaur<-filter(fastSaur, is.element(Masspresent, c(0)))
notumorFoundFastSaur<-notumorFoundFastSaur[notumorFoundFastSaur$proportion_lifespan > 1, ]
tumorFoundFastSaur<-filter(fastSaur, is.element(Masspresent, c(1)))
tumorFoundFastSaur<-tumorFoundFastSaur[tumorFoundFastSaur$proportion_lifespan > 1, ]
ks.test(notumorFoundFastSaur$proportion_lifespan,tumorFoundFastSaur$proportion_lifespan)


fasttumorSaur <- ggplot(fastSaur, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Fast Life History", 
          subtitle="Sauropsids") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

fasttumorSaur

ggsave("fasttumorSaur.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#fast tumor vs no tumor saur

slowSaur <- filter(fastSlowSaur,fastslow == 0)

notumorFoundSlowSaur<-filter(slowSaur, is.element(Masspresent, c(0)))
notumorFoundSlowSaur<-notumorFoundSlowSaur[notumorFoundSlowSaur$proportion_lifespan > 1, ]
tumorFoundSlowSaur<-filter(slowSaur, is.element(Masspresent, c(1)))
tumorFoundSlowSaur<-tumorFoundSlowSaur[tumorFoundSlowSaur$proportion_lifespan > 1, ] 
ks.test(notumorFoundSlowSaur$proportion_lifespan,tumorFoundSlowSaur$proportion_lifespan)


slowtumorSaur<- ggplot(slowSaur, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Masspresent)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Slow Life History",
          subtitle="Sauropsids") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("No Tumor Found", "Tumor Found"))

slowtumorSaur

ggsave("slowtumorSaur.pdf",width=13, height=10, limitsize=FALSE,bg="white")





#tumors found: slow v fast

found<- filter(fastSlow,Masspresent == 1)

TumorFast<-filter(found, is.element(fastslow, c(1)))
TumorFast<-TumorFast[TumorFast$proportion_lifespan > 1, ]
TumorSlow<-filter(found, is.element(fastslow, c(0)))
TumorSlow<-TumorSlow[TumorSlow$proportion_lifespan > 1, ]
ks.test(TumorFast$proportion_lifespan,TumorSlow$proportion_lifespan)

foundSlowMammfast<- ggplot(found, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density Tumors Found",
          subtitle = "All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

foundSlowMammfast

ggsave("foundAll.pdf",width=13, height=10, limitsize=FALSE,bg="white")

#no tumor:slow v fast
nofound<- filter(fastSlow,Masspresent == 0)

noTumorFast<-filter(nofound, is.element(fastslow, c(1)))
noTumorFast<-noTumorFast[noTumorFast$proportion_lifespan > 1, ]
noTumorSlow<-filter(nofound, is.element(fastslow, c(0)))
noTumorSlow<-noTumorSlow[noTumorSlow$proportion_lifespan > 1, ]
ks.test(noTumorFast$proportion_lifespan,noTumorSlow$proportion_lifespan)

noFoundSlowMammfast<- ggplot(nofound, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density No Tumors Found",
          subtitle = "All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

noFoundSlowMammfast

ggsave("nofoundAll.pdf",width=13, height=10, limitsize=FALSE,bg="white")



#tumors found: slow v fast Mammals

foundMamm<- filter(fastSlowMamm,Masspresent == 1)


TumorFastMamm<-filter(foundMamm, is.element(fastslow, c(1)))
TumorFastMamm<-TumorFastMamm[TumorFastMamm$proportion_lifespan > 1, ]
TumorSlow<-filter(foundMamm, is.element(fastslow, c(0)))
TumorSlow<-TumorSlow[TumorSlow$proportion_lifespan > 1, ]
ks.test(TumorFastMamm$proportion_lifespan,TumorSlow$proportion_lifespan)

foundSlowMammfastMamm<- ggplot(foundMamm, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density Tumors Found",
          subtitle = "Mammals") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

foundSlowMammfastMamm

ggsave("foundMamm.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#no tumor:slow v fast Mammals
nofoundMamm<- filter(fastSlowMamm,Masspresent == 0)

noTumorFastMamm<-filter(nofoundMamm, is.element(fastslow, c(1)))
noTumorFastMamm<-noTumorFastMamm[noTumorFastMamm$proportion_lifespan > 1, ]
noTumorSlow<-filter(nofoundMamm, is.element(fastslow, c(0)))
noTumorSlow<-noTumorSlow[noTumorSlow$proportion_lifespan > 1, ]
ks.test(noTumorFastMamm$proportion_lifespan,noTumorSlow$proportion_lifespan)

noFoundSlowMammfastMamm<- ggplot(nofoundMamm, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density No Tumors Found",
          subtitle = "Mammals") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

noFoundSlowMammfastMamm

ggsave("nofoundMamm.pdf",width=13, height=10, limitsize=FALSE,bg="white")





#tumors found: slow v fast Sauropsids

foundSaur<- filter(fastSlowSaur,Masspresent == 1)



TumorFastSaur<-filter(foundSaur, is.element(fastslow, c(1)))
TumorFastSaur<-TumorFastSaur[TumorFastSaur$proportion_lifespan > 1, ]
TumorSlow<-filter(foundSaur, is.element(fastslow, c(0)))
TumorSlow<-TumorSlow[TumorSlow$proportion_lifespan > 1, ]
ks.test(TumorFastSaur$proportion_lifespan,TumorSlow$proportion_lifespan)




foundSlowMammfastSaur<- ggplot(foundSaur, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density Tumors Found",
          subtitle = "Sauropsids") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

foundSlowMammfastSaur

ggsave("foundSaur.pdf",width=13, height=10, limitsize=FALSE,bg="white")


#no tumor:slow v fast Sauropsids
nofoundSaur<- filter(fastSlowSaur,Masspresent == 0)

noTumorFastSaur<-filter(nofoundSaur, is.element(fastslow, c(1)))
noTumorFastSaur<-noTumorFastSaur[noTumorFastSaur$proportion_lifespan > 1, ]
noTumorSlow<-filter(nofoundSaur, is.element(fastslow, c(0)))
noTumorSlow<-noTumorSlow[noTumorSlow$proportion_lifespan > 1, ]
ks.test(noTumorFastSaur$proportion_lifespan,noTumorSlow$proportion_lifespan)

noFoundSlowMammfastSaur<- ggplot(nofoundSaur, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((fastslow)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Fast v. Slow Density No Tumors Found",
          subtitle = "Sauropsids") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c("Slow", "Fast"))

noFoundSlowMammfastSaur


ggsave("nofoundSaur.pdf",width=13, height=10, limitsize=FALSE,bg="white")
