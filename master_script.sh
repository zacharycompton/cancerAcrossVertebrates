#!/bin/bash

# Check required libraries
echo "Checking required R libraries..."
Rscript checkLibrary.R

# Run PGLS analyses
echo "Running PGLS analyses for all life history variables..."
Rscript min20all.R

echo "Running PGLS analyses for mammals..."
Rscript min20mam.R

echo "Running PGLS analyses for aves..."
Rscript min20aves.R

echo "Running PGLS analyses for amphibians..."
Rscript min20amph.R

echo "Running PGLS analyses for sauropsids..."
Rscript min20Sauropsids.R

echo "Running PGLS analyses for vertebrate groups..."
Rscript min20VertGroup.R

# Controlled PGLS analysis
echo "Running controlled PGLS analyses..."
R --slave -e "rmarkdown::render('min20AgeControl.Rmd')"

# Run PGLS for mutation data
echo "Running PGLS for mutation data..."
R --slave -e "rmarkdown::render('mutationRegression.Rmd')"

# Kruskal-Wilcox tests
echo "Running Kruskal-Wilcox tests..."
Rscript kruskalwilcox.R

# Bootstrap simulations
echo "Running bootstrap simulations..."
Rscript bootstrap.R

# Additional PGLS analyses with treatments
echo "Running PGLS for doxirubicin treatment data..."
Rscript doxfunc.R

echo "Running PGLS for carnivore indicators..."
Rscript carnivorePgls.R

# Visualization scripts
echo "Generating visualizations..."
#Rscript superTree.R
#Rscript barTree.R
#Rscript CladeTree.R
#Rscript combined_distributions.R
#Rscript violinsNeoMal.R

# Experimentation scripts
echo "Running experimental scripts..."
#Rscript expExperiments.R

echo "All tasks completed successfully!"
