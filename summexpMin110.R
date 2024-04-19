library(dplyr)

# Read the data from the CSV file
data <- read.csv("415exp.csv")

# Define the desired ranges for the species number
breaks <- c(0, 50, 100, 350)
labels <- c("50", "100", "350")  # Corrected to match the specified ranges

# Add a new column for the species number range
data$speciesnumber_range <- cut(data$n, breaks = breaks, labels = labels, include.lowest = TRUE, right = FALSE)

# Filtering the dataset for rows where min is 1
data_points_min1 <- data[data$min == 1, ]

# Calculating the summary statistics with respect to 'speciesnumber_range'
summary_stats <- data_points_min1 %>%
  group_by(model, speciesnumber_range) %>%
  summarise(
    Total = n(),
    Negative_Slope_Count = sum(pval >= 0.05 & slope < 0, na.rm = TRUE),
    No_Relation_Count = sum(pval <= 0.05 & slope == 0, na.rm = TRUE),
    Positive_Slope_Count = sum(pval >= 0.05 & slope > 0, na.rm = TRUE),
    .groups = 'drop'  # This drops the grouping after summarisation
  ) %>%
  mutate(
    Negative_Slope_Percent = paste0(round((Negative_Slope_Count / Total) * 100, 0), "%"),
    No_Relation_Percent = paste0(round((No_Relation_Count / Total) * 100, 0), "%"),
    Positive_Slope_Percent = paste0(round((Positive_Slope_Count / Total) * 100, 0), "%")
  ) %>%
  select(model, speciesnumber_range, Negative_Slope_Percent, No_Relation_Percent, Positive_Slope_Percent)

# Print the table
print(summary_stats)

write.csv(summary_stats, "./min1experiments.csv")


# Filtering the dataset for rows where min is 1
data_points_min1 <- data[data$min == 10, ]

# Calculating the summary statistics with respect to 'speciesnumber_range'
summary_stats10 <- data_points_min1 %>%
  group_by(model, speciesnumber_range) %>%
  summarise(
    Total = n(),
    Negative_Slope_Count = sum(pval >= 0.05 & slope < 0, na.rm = TRUE),
    No_Relation_Count = sum(pval <= 0.05 & slope == 0, na.rm = TRUE),
    Positive_Slope_Count = sum(pval >= 0.05 & slope > 0, na.rm = TRUE),
    .groups = 'drop'  # This drops the grouping after summarisation
  ) %>%
  mutate(
    Negative_Slope_Percent = paste0(round((Negative_Slope_Count / Total) * 100, 0), "%"),
    No_Relation_Percent = paste0(round((No_Relation_Count / Total) * 100, 0), "%"),
    Positive_Slope_Percent = paste0(round((Positive_Slope_Count / Total) * 100, 0), "%")
  ) %>%
  select(model, speciesnumber_range, Negative_Slope_Percent, No_Relation_Percent, Positive_Slope_Percent)

# Print the table
print(summary_stats10)

write.csv(summary_stats10, "./min10experiments.csv")


