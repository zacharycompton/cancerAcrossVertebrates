# cancerAcrossVertebrates
## Before you start
### Check your library 
```rscript checkLibrary.R```  
You can also run this in RStudio

### checkLibrary.R  
Checks if user has libraries needed to run R scripts

## Data

###  min20-2022.05.16.xlsx 
This is our species level life history data
Contains data and data dictionary
Minimum of 20 records per species

###  min20-2022.05.16.csv
Same life history species data as excel file, just as a csv

### data-description-min20-2022.05.16.txt
Data descrption for min20-2022.05.16.csv

### allrecords.csv
All individual records 

### cleanPath.min20.062822.csv
Individual records with species n >= 20

### fastslow.csv
life history data sorted by slow and fast life history 

### min1.csv
All species level data

### min20DOX.csv
Species level doxirubicin treatment data

### min20RAD_UPDATE.csv
Species level radiation treatment data

### Mutation_Data.csv
Species level mutation rate data

### trophicData.csv
Trophic level data for each species


## Statistical Method Script

### min20all.R
Main PGLS analysis of all life history variables

### min20mam.R
PGLS Analysis of life history variables for mammals

### min20aves.R
PGLS Analysis of life history variables for aves

### min20amph.R
PGLS Analysis of life history variables for amphibians

### min20Sauropsids.R
PGLS Analysis of life history variables for sauropsids

### min20VertGroup.R
PGLS Analysis of life history variables for different vertebrate groups

### min20AgeControl.Rmd
PGLS controlled for available individual age data

### min20GEE.Rmd
Replication of Bull et al.'s phylogenetic life history statistical analysis

### Mltv.AIC.Rmd
Multivariate pgls

### mutationRegression.Rmd
PGLS for mutation data

### pglsPagel.R and pglsPagel.Rmd
Run this code to use pgls in environment. This code is at beginning of all pgls script, so not necessary to run. 

### benignvmal.R
PGLS between benign and malignancy prevalence

### bootstrap.R
Grabs a proportion of population's records and runs PGLS 100 times for both Neoplasia and Malignancy for each significant life history variable.

### carnivorePgls.R
PGLS with Carnivore indicator added to models

### doxfunc.R
PGLS for doxirubicin treatment data

### kruskalwilcox.R
Perform Kruskal and Wilcox test between clades

### arcsinMin20.Rmd || arcsinMin20.html
arcsine tranformation of significant dependent variables in pgls analysis

## Results from Statistical Methods

### bootstrapResult2.0.xlsx
Results of a run of bootstrap.R

### geePglsComp.xlsx
Results of min20GEE.Rmd, min20All.R, and min20AgeControl.Rmd for model comparison

### min20GEEOutput.xlsx
Results of Min20GEE.Rmd

## Script for Visuals

### barTree.R
Creates phylo tree with bars above species name

### CladeTree.R
Ancestral State Reconstruction with Clade bars

### combined_distributions.R
Age of death distribution plots 

### violinsNeoMal.R
Script for violin plots

## Experimentation

### expExperiments.R
Simulation test of pglsPagelSey and compar.GEE. This is a large and and long-running script. 























