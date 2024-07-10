library(dbscan)
library(dplyr)
library(tidyverse)
library(janitor)

# Set the working directory and load the data
csv_file_path <- "C:/Users/ah437/OneDrive/Documents/Work/"
setwd("C:/Users/ah437/OneDrive/Documents/Work")
raw_data_urine <- read_csv("MS Urine.csv")

# Clean column names using janitor
raw_data_urine <- janitor::clean_names(raw_data_urine)

# Optionally, remove specific columns if necessary
raw_data_urine_clean <- raw_data_urine[, -c(3:22)]

# Assuming 'id_and_day' needs to be adjusted by removing the last three characters
raw_data_urine_clean$id_and_day <- substr(raw_data_urine_clean$id_and_day, 1, nchar(raw_data_urine_clean$id_and_day) - 3)

# Remove columns that contain any NA values
raw_data_urine_clean <- raw_data_urine_clean %>% select_if(~!any(is.na(.)))

# names(raw_data_urine_clean)[3:ncol(raw_data_urine_clean)]<-str_split_i(names(raw_data_urine_clean)[3:ncol(raw_data_urine_clean)], "_hmdb", 1)

# Select only numeric columns for scaling
numeric_data <- select_if(raw_data_urine_clean, is.numeric)

# Show data in RStudio's viewer
View(raw_data_urine_clean)
View(numeric_data)

# Calculate mean and standard deviation for each metabolite column grouped by day
stats <- raw_data_urine_clean %>%
  group_by(id_and_day) %>%
  summarise(across(everything(), list(mean = ~mean(.x, na.rm = TRUE), 
                                      sd = ~sd(.x, na.rm = TRUE))))

# View the statistics
stats <- stats %>% select_if(~!any(is.na(.)))

View(stats)

outliers <- raw_data_urine_clean %>%
  group_by(id_and_day) %>%
  mutate(across(where(is.numeric), ~abs(.x - mean(.x, na.rm = TRUE)) > 2 * sd(.x, na.rm = TRUE), .names = "{.col}_outlier"))

# View the dataframe with outliers flagged
View(outliers)

# Remove columns that contain integer data
outliers_cleaned <- outliers %>%
  select(where(~ !is.numeric(.)))

# Adjusting the dataframe to long format, eliminating unnecessary 'type' column
outliers_long <- outliers_cleaned %>%
  pivot_longer(
    cols = -c(id_and_day, day),  # Keep 'id_and_day' and 'day' as identifiers
    names_to = "metabolite",      # Create one new column for metabolite names only
    names_pattern = "(.*)_outlier",  # Regex to extract just the metabolite name before '_outlier'
    values_to = "outlier_status"  # Name of the new column for the TRUE/FALSE values
  )

# View the transformed data to ensure it looks as expected
View(outliers_long)

# Remove rows where 'outlier_status' is FALSE
outliers_long_filtered <- outliers_long %>%
  filter(outlier_status)  # This automatically filters to keep only TRUE values

# View the filtered data to ensure it contains only rows with TRUE in 'outlier_status'
View(outliers_long_filtered)

library(ggplot2)

# Define the directory to save plots
plot_dir <- "C:/Users/ah437/OneDrive/Documents/Work/Metabolite_Plots"
dir.create(plot_dir, showWarnings = FALSE)  # Create the directory if it does not exist, no warning if already exists

# Prepare the list to hold the plots
plot_list <- list()
plot_count <- 0

# Calculate if there are any outliers for each metabolite and group by day
outlier_summary <- outliers %>%
  group_by(id_and_day) %>%
  summarise(across(ends_with("outlier"), ~any(.x == TRUE), .names = "has_{.col}")) %>%
  ungroup()


# Loop through numeric columns and generate plots only for those with outliers
for (metabolite in names(select_if(raw_data_urine_clean, is.numeric))) {
  outlier_column_name <- paste(metabolite, "outlier", sep = "_")
  summary_column_name <- paste("has", outlier_column_name, sep = "_")
  
  # Check if the summary data frame indicates any outliers for this metabolite
  if (any(outlier_summary[[summary_column_name]] == TRUE, na.rm = TRUE)) {
    p <- ggplot(outliers, aes_string(x = "id_and_day", y = metabolite)) +
      geom_line() +
      geom_point(aes_string(color = outlier_column_name)) +
      geom_text(aes(label = day), nudge_y = 0.1, hjust = 0.5, vjust = -1, check_overlap = TRUE, color = "black", size = 3) +
      scale_color_manual(values = c("green", "red")) +
      labs(title = paste("Outliers in", metabolite, "Levels Across Patient ID"),
           x = "ID", y = "Metabolite Level")
    
    # Store the plot in the list for potential further use
    plot_list[[length(plot_list) + 1]] <- p
    
    # Define the filename based on the metabolite and save the plot
    file_name <- paste(plot_dir, paste("Plot_", metabolite, ".png", sep = ""), sep = "/")
    ggsave(file_name, plot = p, width = 10, height = 8, dpi = 300)
  }
}

# Print each plot if it exists
if (length(plot_list) > 0) {
  for (plot in plot_list) {
    print(plot)  # Displays each plot in the RStudio Plots pane
  }
} else {
  print("No outliers found across any metabolites for the days analyzed.")
}

temp<-as_tibble(outliers_long_filtered)
metabolite_count <- temp%>%count(metabolite)%>%arrange(desc(n))
day_count <- temp%>%count(day)
id_count <-  temp%>%count(id_and_day)

write.csv(id_count, paste0(csv_file_path, "outlier_count_id.csv"), row.names = FALSE)
