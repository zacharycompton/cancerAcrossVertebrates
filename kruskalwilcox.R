library(tidyverse)

Data <- read.csv('min20516.csv')

neoClade<-kruskal.test(NeoplasiaPrevalence~Clade,data=Data)


malClade<-kruskal.test(MalignancyPrevalence~Clade,data=Data)


saur <- filter(Data, is.element(Clade, c("Sauropsida")))
mamm <- filter(Data, is.element(Clade, c("Mammalia")))
amph <- filter(Data, is.element(Clade, c("Amphibia")))

#Mammals vs Amphibians neo
mAmphN<-wilcox.test(mamm$NeoplasiaPrevalence, amph$NeoplasiaPrevalence)

#Mammals vs Amphibians mal
mAmphM<-wilcox.test(mamm$MalignancyPrevalence, amph$MalignancyPrevalence)

#Amphibians vs Sauropsids neo
aSaurN<-wilcox.test(amph$NeoplasiaPrevalence, saur$NeoplasiaPrevalence)

#Amphibians vs Sauropsids mal
aSaurM<-wilcox.test(amph$MalignancyPrevalence, saur$MalignancyPrevalence)

mAmphN[3]


PValues<-c(as.numeric(neoClade[3]), as.numeric(malClade[3]),as.numeric(mAmphN[3]),
        as.numeric(mAmphM[3]),as.numeric(aSaurN[3]),as.numeric(aSaurM[3]))
Test<-c("Kruskal Neo",'Kruskal Mal', 'Mammal and Amphbian Neo',"Mammal and Amphibian Mal", 
            "Amphibian and Sauropsid Neo","Amphian and Suaropsid Mal")
results<-data.frame(Test, PValues)
view(results)


