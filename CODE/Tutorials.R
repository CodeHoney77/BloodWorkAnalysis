library(palmerpenguins)
library(ggthemes)
library(tidyverse)
penguins

# Data visualization 

# This data frame contains 8 columns.
# For an alternative view, where you can see all variables and the first few observations of each variable, 
# use glimpse(). Or, if you’re in RStudio, run View(penguins) to open an interactive data viewer.

View(penguins)

# ggplot(data = penguins, mapping = aes(x=flipper_length_mm, y = body_mass_g)) + geom_point() <- scatter plot

ggplot(data = penguins, mapping = aes(x=flipper_length_mm, y = body_mass_g)) + 
  geom_point(mapping= aes(color = species, shape = species)) + 
  geom_smooth(method = "lm") +
  labs(
    title = "Body mass and flipper length",
    subtitle = "Dimensions for Adelie, Chinstrap, and Gentoo Penguins",
    x = "Flipper length (mm)", y = "Body mass (g)",
    color = "Species", shape = "Species"
  ) +
  scale_color_colorblind()

# From this plot we can see that the relationship between body mass and flipper length is positive and linear
 
# Penguin info number of rows in a penguins, number of columns

?penguins #
nrow(penguins)
ncol(penguins)

ggplot(data = penguins, mapping = aes(x=bill_depth_mm , y = bill_length_mm)) + geom_point() 
ggplot(data = penguins, mapping = aes(x=species , y = bill_depth_mm)) + geom_point() 


ggplot(data = penguins, mapping = aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point(mapping = aes(color = bill_depth_mm)) +
  geom_smooth()

# Concise way of plotting
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) + 
  geom_point()

# Visualizing


# Order visualizing bin width matters 

ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram(binwidth = 20)

ggplot(penguins, aes(x = body_mass_g)) +
  geom_histogram(binwidth = 2000) 

# both make it difficult to determine a distribution

ggplot(penguins, aes(x = fct_infreq(species))) +
  geom_bar()

# Geom density makes it easy to determine the distribution of a plot 

ggplot(penguins, aes(y = species)) +
  geom_bar()


# Box plot data and density stuff 

ggplot(penguins, aes(x = species, y = body_mass_g)) +
  geom_boxplot()

ggplot(penguins, aes(x = body_mass_g, color = species)) +
  geom_density(linewidth = 0.75)

ggplot(penguins, aes(x = body_mass_g, color = species, fill = species)) +
  geom_density(alpha = 0.5)

ggplot(penguins, aes(x = island, fill = species)) +
  geom_bar()

ggplot(penguins, aes(x = island, fill = species)) +
  geom_bar(position = "fill")

# 2 variables

ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point()

# 3 Or more variables 
ggplot(penguins, aes(x = flipper_length_mm, y = body_mass_g)) +
  geom_point(aes(color = species, shape = island))

ggplot(data = mpg, mapping = aes(x=hwy, y = displ)) + geom_point(aes(color = model, size = year))

ggplot(penguins, aes(x = bill_depth_mm, y = bill_length_mm)) +
  geom_point(aes(color = species))

# Data transformation

View(flights)

flights |>
  count(arr_delay >120)

longest_delays <- flights |>
  arrange(desc(dep_delay)) 

earliest_flight <- flights |>
  filter(!is.na(dep_time))|>
  arrange(dep_time)

head(earliest_flight)

fastest_flight <- flights |>
  filter(!is.na(air_time), air_time >0) |>
  mutate(speed = distance/ (air_time /60)) |>
  arrange(desc(speed))

library(lubridate)

# Combine year, month, and day into a Date variable and check every day of 2013
all_dates <- flights %>%
  mutate(date = make_date(year, month, day)) %>%
  group_by(date) %>%
  summarise(flights_per_day = n()) %>%
  ungroup()

# Create a sequence of all dates in 2013
full_year_dates <- data.frame(date = seq(from = ymd("2013-01-01"), to = ymd("2013-12-31"), by = "day"))

# Check if there are any days with no flights
missing_days <- anti_join(full_year_dates, all_dates, by = "date")

# Check if there are any missing days
if (nrow(missing_days) > 0) {
  print("There were days with no flights:")
  print(missing_days)
} else {
  print("There was at least one flight every day of 2013.")
}

