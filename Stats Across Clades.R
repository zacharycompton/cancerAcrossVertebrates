
library(dplyr)
library(tidyverse)

# Read the CSV file
data <- read.csv('min20-2022.05.16.csv')

view(data)

# Filter  data for 'Mammalia' in the 'Clade' column
mammalia_data <- filter(data, Clade == 'Mammalia')

# Find  range of 'NeoplasiaPrevalence' values
min_neoplasia_prevalence <- min(mammalia_data$NeoplasiaPrevalence, na.rm = TRUE)
max_neoplasia_prevalence <- max(mammalia_data$NeoplasiaPrevalence, na.rm = TRUE)

# Print range
cat("The range of NeoplasiaPrevalence for Mammalia is from", min_neoplasia_prevalence, "to", max_neoplasia_prevalence)

# Mean for 'NeoplasiaPrevalence' for 'Mammalia'
mean_neoplasia_prevalence <- mean(mammalia_data$NeoplasiaPrevalence, na.rm = TRUE)

# Print mean
cat("The mean of NeoplasiaPrevalence for Mammalia is", mean_neoplasia_prevalence)

# Calculate the median of 'NeoplasiaPrevalence' for 'Mammalia'
median_neoplasia_prevalence <- median(mammalia_data$NeoplasiaPrevalence, na.rm = TRUE)

# Print the median
cat("The median of NeoplasiaPrevalence for Mammalia is", median_neoplasia_prevalence)

# Calculate the mean of 'MalignancyPrevalence' for 'Mammalia'
mean_malignancy_prevalence <- mean(mammalia_data$MalignancyPrevalence, na.rm = TRUE)

# Calculate the median of 'MalignancyPrevalence' for 'Mammalia'
median_malignancy_prevalence <- median(mammalia_data$MalignancyPrevalence, na.rm = TRUE)

# Find the minimum (for range calculation) of 'MalignancyPrevalence' for 'Mammalia'
min_malignancy_prevalence <- min(mammalia_data$MalignancyPrevalence, na.rm = TRUE)

# Find the maximum (for range calculation) of 'MalignancyPrevalence' for 'Mammalia'
max_malignancy_prevalence <- max(mammalia_data$MalignancyPrevalence, na.rm = TRUE)

# Print the results
cat("The mean of MalignancyPrevalence for Mammalia is", mean_malignancy_prevalence, "\n")
cat("The median of MalignancyPrevalence for Mammalia is", median_malignancy_prevalence, "\n")
cat("The range of MalignancyPrevalence for Mammalia is from", min_malignancy_prevalence, "to", max_malignancy_prevalence)


# Filter the data for 'Amphibia' in the 'Clade' column
amphibia_data <- filter(data, Clade == 'Amphibia')

# Calculate the mean of 'NeoplasiaPrevalence' for 'Amphibia'
mean_neoplasia_prevalence_amphibia <- mean(amphibia_data$NeoplasiaPrevalence, na.rm = TRUE)

# Calculate the median of 'NeoplasiaPrevalence' for 'Amphibia'
median_neoplasia_prevalence_amphibia <- median(amphibia_data$NeoplasiaPrevalence, na.rm = TRUE)

# Find the range (minimum and maximum) of 'NeoplasiaPrevalence' for 'Amphibia'
min_neoplasia_prevalence_amphibia <- min(amphibia_data$NeoplasiaPrevalence, na.rm = TRUE)
max_neoplasia_prevalence_amphibia <- max(amphibia_data$NeoplasiaPrevalence, na.rm = TRUE)

# Print the mean, median, and range
cat("The mean of NeoplasiaPrevalence for Amphibia is", mean_neoplasia_prevalence_amphibia, "\n")
cat("The median of NeoplasiaPrevalence for Amphibia is", median_neoplasia_prevalence_amphibia, "\n")
cat("The range of NeoplasiaPrevalence for Amphibia is from", min_neoplasia_prevalence_amphibia, "to", max_neoplasia_prevalence_amphibia)

# Calculate the mean of 'MalignancyPrevalence' for 'Amphibia'
mean_malignancy_prevalence_amphibia <- mean(amphibia_data$MalignancyPrevalence, na.rm = TRUE)

# Calculate the median of 'MalignancyPrevalence' for 'Amphibia'
median_malignancy_prevalence_amphibia <- median(amphibia_data$MalignancyPrevalence, na.rm = TRUE)

# Find the range (minimum and maximum) of 'MalignancyPrevalence' for 'Amphibia'
min_malignancy_prevalence_amphibia <- min(amphibia_data$MalignancyPrevalence, na.rm = TRUE)
max_malignancy_prevalence_amphibia <- max(amphibia_data$MalignancyPrevalence, na.rm = TRUE)

# Print the mean, median, and range
cat("The mean of MalignancyPrevalence for Amphibia is", mean_malignancy_prevalence_amphibia, "\n")
cat("The median of MalignancyPrevalence for Amphibia is", median_malignancy_prevalence_amphibia, "\n")
cat("The range of MalignancyPrevalence for Amphibia is from", min_malignancy_prevalence_amphibia, "to", max_malignancy_prevalence_amphibia)


# Filter the data for 'Sauropsida' in the 'Clade' column
sauropsida_data <- filter(data, Clade == 'Sauropsida')

# Calculate the mean of 'NeoplasiaPrevalence' for 'Sauropsida'
mean_neoplasia_prevalence_sauropsida <- mean(sauropsida_data$NeoplasiaPrevalence, na.rm = TRUE)

# Calculate the median of 'NeoplasiaPrevalence' for 'Sauropsida'
median_neoplasia_prevalence_sauropsida <- median(sauropsida_data$NeoplasiaPrevalence, na.rm = TRUE)

# Find the range (minimum and maximum) of 'NeoplasiaPrevalence' for 'Sauropsida'
min_neoplasia_prevalence_sauropsida <- min(sauropsida_data$NeoplasiaPrevalence, na.rm = TRUE)
max_neoplasia_prevalence_sauropsida <- max(sauropsida_data$NeoplasiaPrevalence, na.rm = TRUE)

# Print the mean, median, and range
cat("The mean of NeoplasiaPrevalence for Sauropsida is", mean_neoplasia_prevalence_sauropsida, "\n")
cat("The median of NeoplasiaPrevalence for Sauropsida is", median_neoplasia_prevalence_sauropsida, "\n")
cat("The range of NeoplasiaPrevalence for Sauropsida is from", min_neoplasia_prevalence_sauropsida, "to", max_neoplasia_prevalence_sauropsida)

# Calculate the mean of 'MalignancyPrevalence' for 'Sauropsida'
mean_malignancy_prevalence_sauropsida <- mean(sauropsida_data$MalignancyPrevalence, na.rm = TRUE)

# Calculate the median of 'MalignancyPrevalence' for 'Sauropsida'
median_malignancy_prevalence_sauropsida <- median(sauropsida_data$MalignancyPrevalence, na.rm = TRUE)

# Find the range (minimum and maximum) of 'MalignancyPrevalence' for 'Sauropsida'
min_malignancy_prevalence_sauropsida <- min(sauropsida_data$MalignancyPrevalence, na.rm = TRUE)
max_malignancy_prevalence_sauropsida <- max(sauropsida_data$MalignancyPrevalence, na.rm = TRUE)

# Print the mean, median, and range
cat("The mean of MalignancyPrevalence for Sauropsida is", mean_malignancy_prevalence_sauropsida, "\n")
cat("The median of MalignancyPrevalence for Sauropsida is", median_malignancy_prevalence_sauropsida, "\n")
cat("The range of MalignancyPrevalence for Sauropsida is from", min_malignancy_prevalence_sauropsida, "to", max_malignancy_prevalence_sauropsida)