setwd("C:/Users/ah437/OneDrive/Documents/Work")

library(tidyverse)
library(ggpubr)
library(rstatix)

# Read and clean the data
raw_data_urine <- read_csv("MS Urine.csv")
raw_data_urine <- janitor::clean_names(raw_data_urine)
raw_data_urine_clean <- raw_data_urine[, -c(3:20)]

raw_data_urine_clean$id_and_day <- substr(raw_data_urine_clean$id_and_day, 1, nchar(raw_data_urine_clean$id_and_day) - 3)

# List of metabolites to analyze
metabolites <- colnames(raw_data_urine_clean)[3:ncol(raw_data_urine_clean)] # Assuming the first two columns are ID and day

# Initialize a data frame to store p-values
results <- data.frame(Metabolite = character(), PValue = numeric(), stringsAsFactors = FALSE)

# Loop through each metabolite
for (metabolite in metabolites) {
  # Perform ANOVA test
  res.aov <- anova_test(data = raw_data_urine_clean, dv = !!sym(metabolite), wid = id_and_day, within = day)
  
  # Extract the ANOVA table
  anova_table <- get_anova_table(res.aov)
  
  # Extract the p-value
  anova_p_value <- anova_table$p[anova_table$Effect == "day"]
  
  # Assign column names based on the ANOVA table, and add metabolite name
  anova_table$Metabolite <- metabolite
  
  # Store the results
  results <- rbind(results, anova_table)
  
  # Plotting
  bxp <- ggboxplot(raw_data_urine_clean, x = "day", y = metabolite, add = "point")
  pwc <- raw_data_urine_clean %>%
    pairwise_t_test(
      formula(paste(metabolite, "~ day")), paired = TRUE,
      p.adjust.method = "bonferroni"
    ) %>%
    add_xy_position(x = "day")
  
  plot <- bxp + 
    stat_pvalue_manual(pwc) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE),
      caption = get_pwc_label(pwc)
    )
  
  print(plot) # Ensure the plot is displayed
}

# Display all collected p-values
View(results)


