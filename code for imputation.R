# Load necessary libraries
library(tidyverse)
library(readr)
library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(raster)
library(data.table)

# Read population data for 2020, 2021, and 2022; 00, 08, 16 for different times
# Find the daytime and nighttime for different countries in UTC
# --- Data Import ---

# Read population data for three consecutive years
pop2020 <- read_csv('C:/.../FB_daily_8hr_pop_tiles_00_2020.csv')
pop2021 <- read_csv('C:/.../FB_daily_8hr_pop_tiles_00_2021.csv')
pop2022 <- read_csv('C:/.../FB_daily_8hr_pop_tiles_00_2022.csv')
pop_combined <- rbind(pop2020, pop2021, pop2022)
remove(pop2020, pop2021, pop2022) # Clear individual year datasets from memory

# Filter and select relevant columns
pop_combined <-  subset(pop_combined, select = -c(density_crisis, density_baseline, clipped_z_score))
pop_combined <- pop_combined %>% filter(!is.na(percent_change))

# Filter rows with 'quadkey' at the same level
quadkey_len <- nchar(pop_combined$quadkey)
length_counts <- table(quadkey_len)
max_quadkey <- names(length_counts[length_counts == max(length_counts)])
pop_quadkey <- pop_combined %>% filter(nchar(pop_combined$quadkey) == max_quadkey)

# Add a 'week' column indicating the day of the week
pop_quadkey <- pop_quadkey %>% mutate(week = lubridate::wday(Date, label = TRUE, abbr = TRUE))   #abbr = FALSE  no abbreviation
pop_baselinenoNA <- pop_quadkey%>% filter(!is.na(n_baseline))

# Calculate the percent of complete records
p_missing <- nrow(pop_baselinenoNA) / nrow(pop_quadkey) 

pop_baseline <- pop_baselinenoNA %>%
  group_by(lat, lon, week) %>%
  summarize(baseline = mean(n_baseline), .groups = 'drop')  # The .groups='drop' prevents the result from being a grouped df


# --- Data Imputation: Step 1A ---
# Substituting missing baseline values with those on other dates
# 1A in the manuscript

# Pivot the data wider for easier imputation
popwide <- pop_baseline %>% 
  pivot_wider(
    names_from = "week",
    values_from = "baseline"
  )

# Calculate mean values of weekdays and weekends for imputation
popwide$weekday <- rowMeans(popwide[c('Mon','Tue','Wed','Thu','Fri')], na.rm=TRUE)
popwide$weekend <- rowMeans(popwide[c('Sun','Sat')], na.rm=TRUE)

# Define workdays
workday <- c('Mon','Tue','Wed','Thu','Fri')

# Replace missing baselines with the mean of workdays(weekends)
# If not all workdays or weekends were missed

for (i in 1:nrow(popwide)) {
  for (j in 1:ncol(popwide)) {
    if (is.na(popwide[i, j])) {
      col_name <- colnames(popwide)[j]
      if (col_name %in% workday) {
        popwide[i, j] <- popwide$weekday[i]
      } else {
        popwide[i, j] <- popwide$weekend[i]
      }
    }
  }
}

# replace missing baseline values of the workdays(weekends) with the mean of weekends (workdays)
# If all workdays(weekends) were missed
popwide[sapply(popwide, is.nan)] <- NA

for (i in 1:nrow(popwide)) {
  for (j in 1:ncol(popwide)) {
    if (is.na(popwide[i, j])) {
      col_name = colnames(popwide)[j]
      if (col_name %in% workday) {
        popwide[i, j] <- popwide$weekend[i]
      } else {
        popwide[i, j] <- popwide$weekday[i]
      }
    }
  }
}

# Prepare the data in long format for merging
popwide_formatch <- subset(popwide, select = -c(weekday, weekend))
popbaselong <- popwide_formatch %>%
  pivot_longer(
    cols = c(3:9),
    names_to = "week",
    values_to = "baseline"
  )

pop_base_1A <- left_join(pop_quadkey, popbaselong, by = c('lat'='lat','lon'='lon',"week" = "week"))

# --- Data Imputation (1A): Step 2 ---
pop_base_1A$num <- pop_base_1A$percent_change/100 *(pop_base_1A$baseline+1) + pop_base_1A$baseline
pop_base_1A$num[pop_base_1A$num < 0] <- 0

# --- Data Imputation: Step 1B ---
# Estimation of missing baseline values for the tiles that lacks all baselines.
# 1B in the manuscript

# Define a function for estimation using worldpop popualtion data
fill_baseline_data <- function(data_a, raster_b, data_c) {
  # Obtain population value from raster image using latitude and longitude
  data_a_sf <- st_as_sf(data_a, coords = c("lon", "lat"), crs = "EPSG:4326")
  raster_data <- raster_b
  filled_data_a <- data_a
  population <- terra::extract(raster_data, data_a_sf)
  filled_data_a$population <- population[,2]
  filled_data_a <- filled_data_a %>% filter(population>0) 
  
  # Find linear relationship between baseline and population size
  models <- setNames(lapply(unique(data_a$week), function(day) {
    a_day <- subset(filled_data_a, week == day)
    lm(baseline ~ population, data = a_day)
  }),unique(data_a$week))
  
  # Obtain population value from raster image using latitude and longitude in dataset c
  data_c_sf <- st_as_sf(data_c, coords = c("lon", "lat"), crs = "EPSG:4326")
  filled_data_c <- data_c
  extracted_c <- terra::extract(raster_data, data_c_sf)
  filled_data_c$population <- extracted_c[,2]
  
  # Use the linear relationship to fill baseline values in dataset c
  for (day in unique(data_c$week)) {
    c_day <- subset(filled_data_c, week == day)
    predictions <- predict(models[[day]], newdata = c_day)
    filled_data_c$est_baseline[data_c$week == day] <- predictions
  }
  # The function processes the input datasets, uses raster data to estimate missing values, and returns the filled data
  return(filled_data_c)
}



# Find missing tiles after imputation 1A
pop_basev <- pop_base_1A %>% filter(!is.na(baseline)) %>% dplyr::select(-quadkey, -n_difference)
pop_basena <- pop_base_1A %>% filter(is.na(baseline)) %>% dplyr::select(-quadkey, -n_difference)

# Calculate the percentage of complete records after step 1A
p_missing1 <- nrow(pop_basev) / nrow(pop_base_1A)

# Load raster data for population in 2019
worldpop_2019 <- rast("C:/.../worldpop/Belgium/bel_ppp_2019_UNadj.tif")

# Resample the WorldPop population data to match the resolution of the FB population data

pop2 <- pop_quadkey 
# find resolution of FB popualtion data
a <- unique(pop2$lon)
b <-unique(pop2$lat)
amin <-min(abs(diff(sort(a))))
bmin <- min(abs(diff(sort(b))))

#raster sample for resample, using 2021-2022 for smaller size
pop_sample <-  pop2 %>%
  group_by(lat,lon) %>%
  summarise(value = mean(percent_change))
min_lat <- min(pop_sample$lat, na.rm = TRUE)
max_lat <- max(pop_sample$lat, na.rm = TRUE)
min_lon <- min(pop_sample$lon, na.rm = TRUE)
max_lon <- max(pop_sample$lon, na.rm = TRUE)

expanded_extent <- terra::ext(min_lon - 0.5, max_lon + 0.5, min_lat - 0.5, max_lat + 0.5)
pop_meanweekrast <- terra::rast(ext=expanded_extent, res=c(amin, bmin), crs='EPSG:4326')


# resample
worldpop_2019 <- aggregate(worldpop_2019, fact=c(amin/worldpop_2019@ptr[["res"]][1], bmin/worldpop_2019@ptr[["res"]][1]), fun=sum, na.rm = TRUE)
worldpopresa_2019 <-  resample(worldpop_2019, pop_meanweekrast , method = "bilinear")
worldpopresa_2019 <- worldpopresa_2019 *(amin/worldpop_2019@ptr[["res"]][1])*(bmin/worldpop_2019@ptr[["res"]][1])

# Use the fill_baseline_data function to estimate missing values in the dataset
data_a <- pop_basev
raster_b <- worldpopresa_2019
data_c <- pop_basena

filled_data <- fill_baseline_data(data_a, raster_b, data_c)

filled_data$baseline <- filled_data$est_baseline
filled_data <- filled_data %>% dplyr::select(-est_baseline,-population)

#set min and max for estimated baseline values (the cut-off Facebook uses for excluding tiles)
filled_data$baseline[filled_data$baseline > 10] <- 10
filled_data$baseline[filled_data$baseline < 0] <- 0

# --- Data Imputation (1B): Step 2 ---
filled_data$num <- filled_data$percent_change/100 *(filled_data$baseline+1) + filled_data$baseline
filled_data$num[filled_data$num < 0] <- 0
filled_data <- filled_data %>% filter(!is.na(num))

# Combine datasets
pop_fin <- rbind(pop_basev, filled_data)

# Calculate complete records percentage after step 1B and 2
p_missing2 <- nrow(pop_fin) / nrow(pop_quadkey)

#correct factor to eliminate the fluctuations in the daily number of observed users
total <- pop_fin %>%
  group_by(Date) %>%
  summarize(num_sum = sum(num)) 

ave <- mean(total$num_sum)
pop_num_est <-  left_join(pop_fin, total, by = c("Date"))
pop_num_est$est <- pop_num_est$num/pop_num_est$num_sum*ave

#export data after imputation
setwd('C:/.../')
fwrite(pop_num_est, "FB_daily_8hr_pop_tiles_00_afterImputation.csv", row.names = F)

