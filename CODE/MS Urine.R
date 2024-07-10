setwd("C:/Users/ah437/OneDrive/Documents/Work")

library(tidyverse)
library(ggpubr)
library(rstatix)

raw_data_urine <- read_csv("MS Urine.csv")
raw_data_urine <-janitor:::clean_names(raw_data_urine)
raw_data_urine_clean <- raw_data_urine[, -c(3:20)]

raw_data_urine_clean$id_and_day <- substr(raw_data_urine_clean$id_and_day, 1, nchar(raw_data_urine_clean$id_and_day) - 3)

View(raw_data_urine_clean)

raw_data_urine_clean %>%
  group_by(day) %>%
  get_summary_stats(creatinine_hmdb0000562, type = "mean_sd") 

raw_data_urine_clean %>%
  group_by(day) %>%
  identify_outliers( creatinine_hmdb0000562)

raw_data_urine_clean %>%
  group_by(day) %>%
  shapiro_test( creatinine_hmdb0000562)


bxp <- ggboxplot(raw_data_urine_clean, x = "day", y="creatinine_hmdb0000562", add="point")
bxp

ggqqplot(raw_data_urine_clean, "creatinine_hmdb0000562", facet.by = "day")

res.aov <- anova_test(data = raw_data_urine_clean, dv =  creatinine_hmdb0000562, wid = id_and_day, within = day)
get_anova_table(res.aov)

# pairwise comparisons
pwc <- raw_data_urine_clean %>%
  pairwise_t_test(
    creatinine_hmdb0000562 ~ day, paired = TRUE,
    p.adjust.method = "bonferroni"
  )

pwc <- pwc %>% add_xy_position(x = "day")

anova_table <- get_anova_table(res.aov)
print(anova_table)

# Extract the p-value directly
anova_p_value <- anova_table$p[anova_table$Effect == "day"]

# Display the p-value
print(anova_p_value)

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "day")
bxp + 
  stat_pvalue_manual(pwc) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
  ) +
  geom_line(aes(group = id_and_day, color = id_and_day), position = position_dodge(0.2))
  geom_point(aes(group = id_and_day, color= id_and_day), position = position_dodge(0.2))
  scale_color_discrete(name= "ID")