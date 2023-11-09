#open libraries
library(gridExtra)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(phytools)
library(cowplot)
library(gridExtra)
library(ggpubr)
library(grid)
library(lemon)
library(sqldf)
#load csv 
Data<-read.csv(file="min20-2022.05.16.csv")
nrow(Data)
oldD<-Data

tree<-read.tree(file="min20Fixed516.nwk")
length(tree$tip.label)

Data$Species <- gsub(" ", "_", Data$Species)
includedSpecies <- Data$Species
pruned.tree<-drop.tip(
  tree, setdiff(
    tree$tip.label, includedSpecies))
length(pruned.tree$tip.label)
pruned.tree <- keep.tip(pruned.tree,pruned.tree$tip.label)
Data$Keep <- Data$Species %in% pruned.tree$tip.label
Data <- Data[!(Data$Keep==FALSE),]


Data <- Data %>%
  mutate(CladeGroup = case_when(
    Clade == "Sauropsida" & Class == "Aves" ~ "Aves",
    Clade == "Sauropsida" & Class == "Reptilia" ~ Orders,
    Clade %in% c("Mammalia", "Amphibia") ~ Clade,
    TRUE ~ NA_character_  # For other cases
  ))


# Assuming your data frame is named 'df'
words_to_keep <- c("Mammalia", "Squamata", "Amphibia", "Aves")

Data <- Data %>%
  filter(grepl(paste(words_to_keep, collapse = "|"), CladeGroup, ignore.case = TRUE))


sum(Data$CladeGroup=="Mammalia")

sum_mammalia_records <- Data %>%
  filter(CladeGroup == "Mammalia") %>%
  summarise(total_mammalia_records = sum(RecordsWithDenominators, na.rm = TRUE))

sum(Data$CladeGroup=="Squamata")

sum_Squamata_records <- Data %>%
  filter(CladeGroup == "Squamata") %>%
  summarise(total_Squamata_records = sum(RecordsWithDenominators, na.rm = TRUE))
sum(Data$CladeGroup=="Amphibia")

sum_Amphibia_records <- Data %>%
  filter(CladeGroup == "Amphibia") %>%
  summarise(total_Amphibia_records = sum(RecordsWithDenominators, na.rm = TRUE))

sum(Data$CladeGroup=="Aves")

sum_Aves_records <- Data %>%
  filter(CladeGroup == "Aves") %>%
  summarise(total_Aves_records = sum(RecordsWithDenominators, na.rm = TRUE))


# Basic violin plot
#Neoplasia
p <- ggplot(Data, aes(x=CladeGroup, y=100*NeoplasiaPrevalence, fill=CladeGroup)) + 
  geom_violin(adjust=1) +
  scale_y_continuous(
    limits = c(0,65))+
  theme(
    legend.text = element_text(size = 10),
  ) +
  scale_fill_manual(values = c("Mammalia" = "#631879FF", "Squamata"= "#008b45ff", "Amphibia"= "#3B4992ff", "Aves" = "red" ),labels=c("98 Species, N=5684", "63 Species, N=2393", "41 Species, N=3061","86 Species, N=4797"))+
  ggtitle("A") +
#  labs(fill = 'Clade')+
  ylab("Neoplasia Prevalence %") +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=0),
  ) +
  theme(plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
print(p)


##Add the jitter points and the mean bar 
pFinal <-p + geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=median, geom="crossbar", size=0.7) + 
  theme_cowplot(12)+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        axis.text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20),
        title=element_text(size = 25),)
print(pFinal)


# Basic violin plot
#Malignancy
m <- ggplot(Data, aes(x=CladeGroup, y=100*MalignancyPrevalence, fill=CladeGroup)) + 
  geom_violin(adjust=1) +
  scale_y_continuous(
    limits = c(0,65))+
  theme(
    legend.text = element_text(size = 10),
  ) +
  scale_fill_manual(values = c("Mammalia" = "#631879FF", "Squamata"= "#008b45ff", "Amphibia"= "#3B4992ff", "Aves" = "red" ),labels=c("98 Species, N=5684", "63 Species, N=2393", "41 Species, N=3061","86 Species, N=4797"))+
  ggtitle("B") +
  #  labs(fill = 'Clade')+
  ylab("Malignancy Prevalence %") +
  theme(
    plot.title = element_text(size = 20, face = "bold")) +
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=0),
  ) +
  theme(plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0))
print(m)


##Add the jitter points and the mean bar 
mFinal <-m + geom_jitter(shape=16, position=position_jitter(0.2))+
  stat_summary(fun=median, geom="crossbar", size=0.7) + 
  theme_cowplot(12)+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size=10),
        axis.text=element_text(size=20),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=20),
        title=element_text(size = 25),)
print(mFinal)


#Arrange both violins into one image
arrange<-ggarrange(pFinal,mFinal,NULL, common.legend = TRUE,legend = "none", heights = c(10,1))

#grab legend function
get_only_legend <- function(plot) {
  plot_table <- ggplot_gtable(ggplot_build(plot))
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box")
  legend <- plot_table$grobs[[legend_plot]]
  return(legend)
}

legend<-get_only_legend(mFinal)

#Reposition legend to middle of page
reposition_legend(arrange, legend = legend, x=.2, y=.05, just = 0)


# get mean and range for neo and mal by clade

#all

alleanNeo<-mean(Data$NeoplasiaPrevalence)
allRangNneo<-range(Data$NeoplasiaPrevalence)

cat("All Neoplasia Mean: ", alleanNeo, "All Neoplasia Range", allRangNneo)

alleanMal<- mean(Data$MalignancyPrevalence)
mamRangeMal<-range(Data$MalignancyPrevalence)

cat("All Malignancy Mean: ", alleanMal, "All Malginancy Range", mamRangeMal)


#Mammals
Mamm<- filter(Data, is.element(Clade, c("Mammalia")))
mamMeanNeo<-mean(Mamm$NeoplasiaPrevalence)
mammRangNneo<-range(Mamm$NeoplasiaPrevalence)

cat("Mammmal Neoplasia Mean: ", mamMeanNeo, " Mammal Neoplasia Range", mammRangNneo)

mamMeanMal<- mean(Mamm$MalignancyPrevalence)
mamRangeMal<-range(Mamm$MalignancyPrevalence)

cat("Mammmal Malignancy Mean: ", mamMeanMal, " Mammal Malginancy Range", mamRangeMal)

#Sauropsids
Saur<- filter(Data, is.element(Clade, c("Sauropsida")))
saurMeanNeo<-mean(Saur$NeoplasiaPrevalence)
saurRangNneo<-range(Saur$NeoplasiaPrevalence)

cat("Sauropsid Neoplasia Mean: ", saurMeanNeo, " Sauropsid Neoplasia Range", saurRangNneo)

saurMeanMal<- mean(Saur$MalignancyPrevalence)
saurRangeMal<-range(Saur$MalignancyPrevalence)

cat("Sauropsid Malignancy Mean: ", saurMeanMal, " Sauropsid Malginancy Range", saurRangeMal)

#Amphibia
amph<- filter(Data, is.element(Clade, c("Amphibia")))
amphMeanNeo<-mean(amph$NeoplasiaPrevalence)
amphRangNneo<-range(amph$NeoplasiaPrevalence)

cat("Amphibian Neoplasia Mean: ", amphMeanNeo, " Sauropsid Neoplasia Range", amphRangNneo)

amphMeanMal<- mean(amph$MalignancyPrevalence)
amphRangeMal<-range(amph$MalignancyPrevalence)

cat("Amphibian Malignancy Mean: ", amphMeanMal, " Sauropsid Malginancy Range", amphRangeMal)








