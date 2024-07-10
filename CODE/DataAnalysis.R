library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(stats)
library(stringr)

# Defining variables, you will have to change file path depending on where you saved them 
MSBlood <- read_excel("C:/Users/ah437/OneDrive/Documents/Work/MasterDOMINO.xlsx", sheet = "MS Blood")
MSUrine <- read_excel("C:/Users/ah437/OneDrive/Documents/Work/MasterDOMINO.xlsx", sheet = "MS Urine")

# Clean the dataset 
names(MSUrine)[1] <- "ID_and_DAY"
names(MSUrine)[2] <- "Day"
names(MSBlood)[1] <- "ID_and_Day"
names(MSBlood)[3] <- "Day"
names(MSBlood)[4] <- "ID"
names(MSBlood)[23] <- "Order_of_sample_run"
names(MSBlood)[20] <- "severityScore"
names(MSBlood)[21] <- "mentalFatigue"

# Use the View command to see any dataset 
# View(data)
# rm(data) this will delete data can be used to reload data if you altered it. 

# General patient data
patientData <- MSBlood[, c(1:22)]
View(MSBlood)

ggplot(data = MSBlood, mapping = aes(x= `AcylCarnitine 12:0`, y = `AcylCarnitine 13:0`)) + 
  geom_point(mapping = aes(color = Day, shape = `(Control, treated etc.)`))

filtered_data <- MSBlood %>%
  filter(Day %in% c(1, 4))

All_delta_changes_D4_D1 <- filtered_data %>%
  group_by(ID, `(Control, treated etc.)`) %>%
  arrange(ID, Day) %>%
  summarise(across(.cols = 22:(ncol(.)-2), ~.[Day == 4] - .[Day == 1], .names = "Delta_{.col}"), .groups = "drop")

View(All_delta_changes_D4_D1)

average_deltas <- All_delta_changes_D4_D1 %>%
  group_by(`(Control, treated etc.)`) %>%
  summarise(across(starts_with("Delta"), mean, na.rm = TRUE))

View(average_deltas)

difference_results <- average_deltas %>%
  summarise(across(starts_with("Delta"), ~ .[1] - .[2], .names = "Difference_{.col}"))

# Print the results

View(difference_results)

# Modify the ID_and_DAY column to remove the day part
MSUrine$ID_and_DAY <- str_extract(MSUrine$ID_and_DAY, "^[^-]+-[^-]+")

# Define day pairs for comparison
day_pairs <- list(
  c("D4", "D2"), 
  c("D5", "D4"), 
  c("D8", "D7"), 
  c("D6", "D5"),
  c("D7", "D6"), 
  c("D8", "D2")
)

# Function to calculate differences for given days
calculate_differences <- function(day2, day1) {
  # Filter data for the specified days
  data_day1 <- MSUrine %>% filter(Day == day1) %>% select(ID_and_DAY, everything()[21:ncol(MSUrine)])
  data_day2 <- MSUrine %>% filter(Day == day2) %>% select(ID_and_DAY, everything()[21:ncol(MSUrine)])
  
  # Ensure both dataframes are in the same order
  data_day1 <- data_day1[order(data_day1$ID_and_DAY),]
  data_day2 <- data_day2[order(data_day2$ID_and_DAY),]
  
  # Calculate changes by subtracting the values of the first day from the second day
  changes <- data_day2[-1] - data_day1[-1]  # Excluding the ID column from subtraction
  
  # Combine the changes with the patient IDs
  results <- data.frame(ID_and_DAY = data_day1$ID_and_DAY, changes)
  names(results)[1] <- "ID"
  
  return(results)
}

# Apply the function to each pair of days and store results
results_list <- lapply(day_pairs, function(pair) calculate_differences(pair[1], pair[2]))

# Optionally, combine all results into one dataframe with an additional column specifying the day comparison
all_results <- do.call(rbind, lapply(seq_along(results_list), function(i) {
  day_comparison <- paste(day_pairs[[i]][1], day_pairs[[i]][2], sep="-")
  cbind(ID = results_list[[i]]$ID, Day_Comparison = day_comparison, results_list[[i]][-1])
}))

# View the final results
View(all_results)