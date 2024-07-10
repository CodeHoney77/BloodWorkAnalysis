library(dbscan)
library(dplyr)
library(tidyverse)
library(janitor)

# Set the working directory and load the data
setwd("C:/Users/ah437/OneDrive/Documents/Work")
raw_data_urine <- read_csv("MS Urine.csv")
raw_data_urine <- janitor::clean_names(raw_data_urine)
raw_data_urine_clean <- raw_data_urine[, -c(3:20)]

# Select only numeric columns for scaling
numeric_data <- select_if(raw_data_urine_clean, is.numeric)

# Apply scaling
metabolite_data_scaled <- scale(numeric_data)

# Assuming 'id_and_day' is the column that identifies each sample
identifiers <- raw_data_urine_clean['id_and_day']

# Now proceed with selecting and scaling the numeric data
numeric_data <- select_if(raw_data_urine_clean, is.numeric)
metabolite_data_scaled <- scale(numeric_data)

# Convert scaled data back to a data frame and reattach the identifier
metabolite_data_scaled <- as.data.frame(metabolite_data_scaled)
metabolite_data_scaled <- bind_cols(identifiers, metabolite_data_scaled)

# Calculate LOF scores and add them to the data frame
k_value <- ceiling(log(nrow(metabolite_data_scaled)))
lof_scores <- lof(metabolite_data_scaled[, -1], k = k_value)  # Exclude identifier column from LOF calculation
metabolite_data_scaled$LOF_Score = lof_scores

# View the data to check identifiers and LOF scores
View(metabolite_data_scaled)

# Define a threshold for outliers
threshold = 1.5  # Adjust this based on your dataset and analysis needs

# Filter outliers based on the LOF score
outliers = filter(metabolite_data_scaled, LOF_Score > threshold)

# View outliers, now including the identifiers
View(outliers)

metabolite_data_scaled$outlier <- ifelse(lof_scores > threshold, "Outlier", "Inlier")

View(metabolite_data_scaled)

file_path <- "C:/Users/ah437/OneDrive/Documents/Work/metabolite_data_scaled.csv"
write.csv(metabolite_data_scaled, file_path, row.names = FALSE)
