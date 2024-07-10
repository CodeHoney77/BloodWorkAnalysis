setwd("C:/Users/ah437/OneDrive/Documents/Work")

library(tidyverse)
library(janitor)
library(ggpubr)
library(rstatix)
library(factoextra)
library(pheatmap)
library(qcc)
library(MASS)
library(dplyr)

# Load and clean the data
raw_data_urine <- read_csv("MS Urine.csv")
raw_data_urine <- janitor::clean_names(raw_data_urine)
raw_data_urine_clean <- raw_data_urine[, -c(3:20)]
raw_data_urine_clean$id_and_day <- substr(raw_data_urine_clean$id_and_day, 1, nchar(raw_data_urine_clean$id_and_day) - 3)

View(raw_data_urine_clean)

# Exclude the first two columns and keep all others
metabolite_data <- raw_data_urine_clean[, sapply(raw_data_urine_clean, is.numeric)]
standard_deviations <- apply(raw_data_urine_clean[, sapply(raw_data_urine_clean, is.numeric)], 2, sd ,na.rm= TRUE)
means <- apply(raw_data_urine_clean[, sapply(raw_data_urine_clean, is.numeric)], 2, mean, na.rm = TRUE)

View(means)
View(standard_deviations)
View(metabolite_data)

# Convert means and standard deviations into a dataframe
results_df <- data.frame(
  Metabolite = names(means),
  Mean = means,
  Standard_Deviation = standard_deviations
)

View(results_df)

# Calculate lower and upper bounds
lower_bounds <- sweep(metabolite_data, 2, means - 2 * standard_deviations, FUN = "-")
upper_bounds <- sweep(metabolite_data, 2, means + 2 * standard_deviations, FUN = "+")

outlier_matrix <- metabolite_data < lower_bounds | metabolite_data > upper_bounds

# Initialize a list to hold data frames for each metabolite's outliers
outlier_values <- lapply(seq_along(metabolite_data), function(i) {
  # Extract indices where the outlier condition is true
  outlier_indices <- which(outlier_matrix[, i])
  if (length(outlier_indices) > 0) {
    # Construct a data frame for outliers in this metabolite
    data.frame(
      ID = outlier_indices,
      Metabolite_Value = metabolite_data[outlier_indices, i],
      Metabolite_Name = names(metabolite_data)[i]
    )
  } else {
    NULL  # Return NULL if no outliers, to avoid errors in combining data frames
  }
})

# Combine all outlier data frames into a single data frame
all_outliers <- do.call(rbind, outlier_values)
View(all_outliers)

# Assuming 'raw_data_urine_clean' is your dataset already loaded
dayDifference <- raw_data_urine_clean %>%
  group_by(id_and_day) %>%
  mutate(
    previous_day = lag(day),  # Capture the previous day
    across(3:ncol(.)-1, ~ . - lag(.), .names = "difference_{.col}")  # Calculate differences for all metabolite columns
  ) %>%
  filter(!is.na(previous_day)) %>%  # Filter out rows where there is no previous day (first entry of each group)
  mutate(
    day_pair = ifelse(!is.na(previous_day), paste(day, previous_day, sep = "-"), NA)  # Create day pair notation
  ) %>%
  dplyr::select(id_and_day, day_pair, starts_with("difference_"))  # Select the id, day_pair, and all difference columns

# Output results
View(dayDifference)

# Assuming 'raw_data_urine_clean' is your dataset already loaded
metabolite_variances <- raw_data_urine_clean %>%
  dplyr::select(3:ncol(.)) %>%  # Select only metabolite columns, starting from the third column
  summarise(across(everything(), ~var(.x, na.rm = TRUE)))  # Calculate variance for each column, removing NA values using the updated syntax

# Output results
View(metabolite_variances)

# Scale and Center the Metabolite Data
dataScaled <- raw_data_urine_clean %>% 
  dplyr::select(-id_and_day, -day) %>%
  scale(center= TRUE, scale = TRUE)

View(dataScaled)

pca_result <- prcomp(dataScaled, center = TRUE, scale. = TRUE)

pca_summary <- summary(pca_result)
var_explained <- pca_summary$importance[2,] * 100  # Convert to percentage

pca_data <- as.data.frame(pca_result$x) %>%
  mutate(Individual = raw_data_urine_clean$id_and_day, Day = raw_data_urine_clean$day)

# Plotting by Individual with variance explained in the title
ggplot(pca_data, aes(x = PC1, y = PC2, color = Individual)) +
  geom_point(aes(shape = Individual), size = 3) +
  theme_minimal() +
  ggtitle(sprintf("PCA of Metabolites by Individual (PC1: %.2f%%, PC2: %.2f%%)",
                  var_explained['PC1'], var_explained['PC2']))

# Plotting by Day with variance explained in the title
ggplot(pca_data, aes(x = PC1, y = PC2, color = Day)) +
  geom_point(aes(shape = Day), size = 3) +
  theme_minimal() +
  ggtitle(sprintf("PCA of Metabolites by Day (PC1: %.2f%%, PC2: %.2f%%)",
                  var_explained['PC1'], var_explained['PC2']))

# For Individual
fviz_pca_ind(pca_result, label = "none", habillage = pca_data$Individual,
             addEllipses = TRUE, ellipse.level = 0.95) +
  ggtitle("PCA of Metabolites by Individual with Ellipses")

# For Day
fviz_pca_ind(pca_result, label = "none", habillage = pca_data$Day,
             addEllipses = TRUE, ellipse.level = 0.95) +
  ggtitle("PCA of Metabolites by Day with Ellipses")