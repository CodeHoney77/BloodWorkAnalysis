library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)

# Defining variables, you will have to change file path depending on where you saved them 
MSBlood <- read_excel("C:/Users/ah437/OneDrive/Documents/Work/MasterDOMINO.xlsx", sheet = "MS Blood")

# Clean the dataset 
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
  filter(Day %in% c(1, 4, 14))

All_delta_changes_D4_D1 <- filtered_data %>%
  group_by(ID, `(Control, treated etc.)`) %>%
  arrange(ID, Day) %>%
  summarise(across(.cols = 23:(ncol(.)-2), ~.[Day == 4] - .[Day == 1], .names = "Delta_{.col}"), .groups = "drop")

View(All_delta_changes_D4_D1)

average_deltas <- All_delta_changes_D4_D1 %>%
  group_by(`(Control, treated etc.)`) %>%
  summarise(across(starts_with("Delta"), mean, na.rm = TRUE))

View(average_deltas)

difference_results <- average_deltas %>%
  summarise(across(starts_with("Delta"), ~ .[1] - .[2], .names = "Difference_{.col}"))

# Print the results
print(difference_results)

View(difference_results)
