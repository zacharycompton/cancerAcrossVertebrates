# Cancer Across Vertebrates

## min20all.R
This file runs PGLS (Phylogenetic Generalized Least Squares) models for life history traits across all species. When running the code in R, you must have all packages installed.
Regression plots are also outputted.

Use this line of code to install a missing package: 

```
install.packages("phytools")
```

This will need to be done for packages you may not have installed


## min20mam.R
This file also preforms PGLS, just for mammals and their available life history traits

## barTree.R
This file creates the large bar graph tree you see in the Cancer Across Vertebrates paper

## CladeTree.R
This file creates the Ancestral State Reconstructed Phylogeny. This file outputs trees for each clade, then creates labels for them.

## combined_distributions.R
This file creates the density age plots.

## fastvslowdensity.R 
This file creates density plots for multiple groupings of fast and slow life history species.

## doxfunc.R 
This file creates functional PGLS regressions and plots for species under a doxorubicin regimen.

## radfunc.R
This file creates functional PGLS regressions and plots for species under a radiation regimen.

## mutationRegression.R
This creates a PGLS regressino and plot using mutation data

## Mltv.AIC.Rmd

This creates pgls multivariate regressions for multiple life history traits.

## kruskalwilcox.R
This file carries out kruskal and wilcox tests between and across clades


