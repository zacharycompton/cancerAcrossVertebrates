#open libraries
library(patchwork)
library(tidyverse)
library(cowplot)
library(ggrepel)


#read csv
cleanPath <- read.csv("cleanPath.min20.062822.csv")
min20Data <- read.csv("min20516.csv")
view(min20Data)

#linear model mav vs benign
malben<-lm(MalignancyPrevalence~BenignPrevalence,data = min20Data) 
summary(malben)
summary(malben)$r.squared
summary(malben)$coefficients[2,4]

#ks test for benign lifespan vs malignant lifespan
benign<-filter(cleanPath, is.element(Malignant, c(0)))
malignant<-filter(cleanPath, is.element(Malignant, c(1)))
ks.test(benign$proportion_lifespan,malignant$proportion_lifespan)
cleanPath<- filter(cleanPath,Malignant >= 0)


#grab r squared
r2<-signif(as.numeric(summary(malben)$r.squared)^2, digits= 2)
pval<-signif(as.numeric(summary(malben)$coefficients[2,4])^2, digits= 2)



#plot malignancy prev vs benign prev
ggplot(min20Data, aes(y=MalignancyPrevalence*100, x= BenignPrevalence*100))+
  scale_color_manual(values = c("Mammalia" = "#631879FF", "Sauropsida"= "#008b45ff", "Amphibia"= "#3B4992ff" ),)+
  scale_y_continuous(
    limits = c(0,50),
    breaks = c(0, 25,50),
    labels = c(0, 25,50))+
  scale_x_continuous(
    limits = c(0,30),
    breaks = c(0, 15,30),
    labels = c(0, 15,30)
  )+
  geom_abline(intercept = coef(malben)[1], slope =  coef(malben)[2],
              color = 'grey',size = 1.2) +
  theme_cowplot(12)+
  theme(axis.title = element_text(size = 18))+
  ylab("Malignancy Prevalence (%)") +
  xlab("Benign Prevalence (%)") +
  geom_point(aes(colour= Clade, size = RecordsWithDenominators)) +
  scale_size(name   = "Total Necropsies",
             breaks = c(20,100,200,300,477),
             labels =  c(20,100,200,300,477))+
  geom_text_repel(aes(label=ifelse( MalignancyPrevalence > .3,as.character(common_name),'')),max.overlaps = Inf,size=5, direction = "y")+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  labs(title = "Malignancy Prevalence v. Benign Prevalence",
       subtitle =bquote(p-value:.(pval)~R^2:.(r2)))+
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(legend.position = "bottom")+
  labs(colour="Clade", size="Total Necropsies")



#save 
ggsave(filename='benignvmal.pdf', width=13, height=10, limitsize=FALSE,bg="white")



#plot prop lifespan malignant factor
bvm<- ggplot(cleanPath, aes(x=proportion_lifespan*100, y=..scaled..,fill=factor((Malignant)))) + 
  geom_density(alpha=0.25) + 
  ggtitle("Tumor v. No Tumor Slow Life History",
          subtitle="All Species") +
  xlab("Age at Death as a Percentage of Species Lifespan") + ylab("Normalized Frequency") + 
  coord_cartesian(xlim = c(0,150),ylim = c(0,1))+theme_cowplot(12)+ theme(legend.position = "bottom")+
  scale_fill_discrete(name="",labels=c( "No Tumor Found","Benign", "Malignant"))

bvm