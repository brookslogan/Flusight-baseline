
##The following code produces the peak week and peak intensity historic baseline estimates for the 2024/2025 season of the FluSight Challenge.

##Packages
library(lubridate)
library(readr)
library(dplyr)
library(tidyr)
library(checkmate)

#specify output path
userid <- Sys.info()["user"]
output_dirpath <- paste0("C:/Users/",userid,"/Desktop/GitHub/Flusight-baseline/weekly-submission/forecasts/Flusight-seasonal-baseline/")
hub_outputpath <- paste0("C:/Users/",userid,"/Desktop/Github/FluSight-forecast-hub/model-output/Flusight-seasonal-baseline/")

#########################
## Functions
# Reference Date
curr_else_next_date_with_ltwday <- function(date, ltwday) {
  assert_class(date, "Date")
  assert_integerish(ltwday, lower = 0L, upper = 7L)
  #
  date + (ltwday - as.POSIXlt(date)$wday) %% 7L
}
forecast_as_of_date <- Sys.Date()
reference_date <- curr_else_next_date_with_ltwday(forecast_as_of_date, 6L) # Saturday

# Conversion from epi week to date
epiweekToDate <- function(year, weekno) {
  # Create a date object for January 1 of the given year
  first_day_of_year <- as.Date(paste(year, "-01-01", sep=""))
  
  # Calculate the first Monday of the given ISO week (week starts on Monday)
  start_date <- first_day_of_year + weeks(weekno - 1)
  
  # Adjust the start_date to the first Monday of the given week
  weekday_start <- wday(start_date, week_start = 1)  # Week starts on Monday
  
  # If the start date is not Monday, calculate the first Monday of the week
  if (weekday_start != 1) {
    start_date <- start_date - days(weekday_start - 1)
  }
  
  # Now that we know the first Monday, adjust to Saturday (5 days after Monday)
  last_saturday <- start_date + days(5)  # Saturday is 5 days after Monday
  
  return(last_saturday)
}


#########################
## Load Data

#FluSurv-NET data 2010/11 through 2019/20
flusurv_state <- read.csv(paste0("C:/Users/",userid,"/Desktop/Github/Flusight-baseline/seasonal-historic/FluSurv-NET_pastdata.csv"))
fsn_data <- flusurv_state %>% 
  filter(AGE.CATEGORY=="Overall", SEX.CATEGORY=="Overall", RACE.CATEGORY=="Overall", VIRUS.TYPE.CATEGORY=="Overall") %>% 
  mutate(weekly_rate = as.numeric(WEEKLY.RATE)) %>% 
  rename("season"="YEAR", "epidemic_week"="WEEK", "location"= "CATCHMENT") %>% 
  filter(epidemic_week <= 22 | epidemic_week >= 35) %>% 
  mutate(`epidemic_week`=as.numeric(`epidemic_week`), week = ifelse(epidemic_week < 36, epidemic_week + 52, epidemic_week))%>% 
  filter(!(season %in% c("2009-10","2020-21","2021-22", "2023-24", "2024-25"))) %>% 
  select(season, epidemic_week, week, weekly_rate, location)


# State population data
state_pop<-read.csv(paste0("C:/Users/",userid,"/Desktop/GitHub/FluSight-forecast-hub/auxiliary-data/locations.csv")) %>% 
  select(location, abbreviation, population)
# list of locations codes
states <- state_pop %>% pull(location)


# NHSN data 2021-2024
nhsn <- read_csv("../Documents/Data/Flu Admissions 21-24.csv")

# Format NHSN data
nhsn_rates <- nhsn %>%
  mutate(
    date = mdy(date),  # Convert MM/DD/YYYY to Date
    epidemic_week = epiweek(date),  # Extract epidemic week using lubridate
    year = year(date)) %>% 
    filter(epidemic_week <= 22 | epidemic_week >= 35) %>% 
    mutate(week = ifelse(epidemic_week < 36, epidemic_week + 52, epidemic_week)) %>% 
  group_by(state, year, epidemic_week,week) %>% # Aggregate by state and week
  summarise(
    total_flu_admissions = sum(admissions, na.rm = TRUE),
    .groups = 'drop'  # To drop the grouping after summarization
  )%>%
  mutate(
    season = case_when(
      epidemic_week >= 35 & epidemic_week <= 52 ~ paste(year, year + 1, sep = "/"),  # For weeks 35-52 of a year
      epidemic_week >= 1 & epidemic_week <= 20 ~ paste(year - 1, year, sep = "/"),  # For weeks 1-20 of the following year
      TRUE ~ NA_character_  # For weeks outside of the flu season (unlikely, but safe to include)
    )
  ) %>% filter(!is.na(season)) %>% filter(!(season %in% c("2020/2021","2019/2020"))) %>% 
  left_join(state_pop, by=c("state"="abbreviation")) %>%  # Calculate rates using state populations
  mutate(rate=total_flu_admissions/population *100000) %>% select(epidemic_week, week, season, rate, state) %>% 
  rename("weekly_rate"= "rate", "location"="state")

#Combine FSN and NHSN data
combo<- rbind(fsn_data, nhsn_rates) %>% filter(!is.na(weekly_rate))



####################
##Peak Week Function
predict_peak_week <- function(these_data) {
  
  # Create data table of peak weeks across years
  these_peaks <- these_data %>%
    group_by(season, location) %>%
    mutate(observation = round(weekly_rate, 1)) %>%
    filter(observation == max(observation, na.rm = TRUE)) %>%
    ungroup() 
  
  # Generate function approximating kernel of past peak values
  pkwk_kernel <- these_peaks %>%
    pull(week) %>%
    density(kernel = "gaussian", bw = "sj", from = min(these_data$week)) %>%
    approxfun(rule = 1:2)
  
  pkper_kernel <- these_peaks %>%
    pull(weekly_rate) %>%
    density(kernel = "gaussian", bw = "sj", from = min(these_data$weekly_rate)) %>%
    approxfun(rule = 1:2)
  
  #### Peak week forecasts
  pkwk_pred <- tibble(reference_date = rep(reference_date, 35),
                      target = rep("peak week inc flu hosp", 35),
                      horizon = rep(NA,35),
                      target_end_date = rep(NA, 35),
                      output_type = rep("pmf", 35),
                      unit = rep("week", 35),
                      bin_start_incl = seq(40, 74, 1),
                      bin_end_notincl = seq(41, 75, 1))%>%
    rowwise() %>%
    mutate(value = integrate(pkwk_kernel, bin_start_incl, bin_end_notincl)[[1]]) %>%
    ungroup() %>%
    mutate(bin_start_incl = ifelse(bin_start_incl > 52, bin_start_incl - 52,
                                   bin_start_incl),
           bin_end_notincl = ifelse(bin_end_notincl > 53, bin_end_notincl - 52,
                                    bin_end_notincl),
           year = ifelse(bin_start_incl >= 40, 2024, 2025),  # Assign year based on epidemic week
           output_type_id = as.Date(mapply(epiweekToDate, year, bin_start_incl)),
           value = value / sum(value)) %>% 
    rename(epidemic_week = bin_start_incl) %>% 
    select(-c(bin_end_notincl, year, unit)) %>% 
    select(-value, everything(), value)
    
  return(pkwk_pred)
}

peak_week <- predict_peak_week(combo)



#########################
# Redistributing Probabilities
  ##The distribution was calculated across weeks 40-52 and 1-22 to encompass a typical influenza season, and the below code
  ##redistributes the probabilities of weeks not included in the current season's FluSight Challenge to the remaining weeks.

  # Define the weeks to remove (weeks 40-46)
  removed_weeks <- 40:46
  # Define the remaining weeks (weeks 47-52 and 1-22)
  remaining_weeks <- c(47:52, 1:22)

#Function for reallocation
redistribute_probabilities <- function(data, removed_weeks, remaining_weeks) {
  
  # Step 1: Calculate the total probability of the removed weeks (40-46)
  removed_prob_sum <- sum(data$value[data$epidemic_week %in% removed_weeks])
  
  # Step 2: Calculate the total probability of the remaining weeks (47-52, 1-22)
  remaining_prob_sum <- sum(data$value[data$epidemic_week%in% remaining_weeks])
  
  # Step 3: Proportional redistribution of the removed probability across the remaining weeks
  # Calculate the proportional redistribution for each remaining week
  adjusted_remaining_probs <- data$value[data$epidemic_week %in% remaining_weeks] + 
    (removed_prob_sum * data$value[data$epidemic_week %in% remaining_weeks] / remaining_prob_sum)
  
  # Step 4: Update the probabilities for the remaining weeks
  data$value[data$epidemic_week %in% remaining_weeks] <- adjusted_remaining_probs
  
  # Step 5: Remove remaining weeks and epidemic week column
  data <- data %>% filter(epidemic_week!= removed_weeks) %>% 
    select(-epidemic_week)
  
  # Return updated data
  return(data)
}

peak_week <- redistribute_probabilities(peak_week, removed_weeks, remaining_weeks)

# Repeat for each state
peak_week <- peak_week %>% 
  crossing(location = states) %>% arrange(location) %>% 
  mutate(output_type_id = as.character(output_type_id))



#########################
## Peak Intensity

#specify quantiles 
quantiles <- c(0.01, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 
               0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 
               0.95, 0.975, 0.99)

calculate_kde_quantiles <- function(data){
  
# Loop through each week, but skip those with fewer than 2 data points
kde_results <- list()

for (week_num in 40:72) {
  # Filter data for the current week
  week_data <- data %>% filter(week == week_num)
  
  # Skip weeks with fewer than 2 data points
  if (nrow(week_data) < 2) {
    next  # Skip this week if not enough data points
  }
  
  # Ensure that WEEKLY.RATE is numeric
  week_data$weekly_rate <- as.numeric(week_data$weekly_rate)
  
  # Fit the kernel density estimate (bandwidth chosen automatically)
  kde <- tryCatch(
    density(week_data$weekly_rate, bw = "SJ"),
    error = function(e) NULL  # Catch errors and assign NULL if the calculation fails
  )
  
  # If the KDE was successful, store the result
  if (!is.null(kde)) {
    kde_results[[week_num]] <- kde
  }
}

# Create a data frame to store the quantile results
quantile_predictions <- data.frame(
  week = rep(40:72, each = length(quantiles)),
  quantile = rep(quantiles, times = length(40:72)),
  value = NA
)

# Loop through each filtered week and calculate the quantiles from the KDE
for (week_num in 40:72) {
  # Get the KDE for the current week
  kde <- kde_results[[week_num]]
  
  # Check if there is a valid KDE
  if (!is.null(kde)) {
    # Calculate the quantiles based on the KDE's cumulative distribution
    for (q in quantiles) {
      # Calculate the quantile from the cumulative distribution
      quantile_predictions$value[quantile_predictions$week == week_num & 
                                      quantile_predictions$quantile == q] <- 
        quantile(kde$x, probs = q)
    }
  }
}

# Calculate the peak intensity for each quantile
peak_intensity <- quantile_predictions %>%
  group_by(quantile) %>%
  summarize(value = max(value, na.rm = TRUE)) %>% # Handle NA values
  rename("output_type_id"="quantile") %>% 
  mutate(reference_date=reference_date,
         target = "peak inc flu hosp", 
         horizon = NA,
         target_end_date = NA,
         output_type = "quantile") %>%
  # Expand the grid by repeating the original data for each state
  crossing(location = states) %>% 
  arrange(location)%>% 
  left_join(state_pop, by = c("location" = "location")) %>% 
  mutate(value = round((value * population) / 100000),
         output_type_id = as.character(output_type_id)) %>% 
  select(-c(population, abbreviation)) %>% 
  select(reference_date, target,horizon, target_end_date, location, output_type, output_type_id, value)

# Return peak intensity
return(peak_intensity = peak_intensity)
}

peak_intensity <- calculate_kde_quantiles(combo)



##Export csv

full <- bind_rows(peak_intensity, peak_week) %>% mutate(location=as.character(location))

write.csv(full, file = paste0(output_dirpath,sprintf("%s-FluSight-seasonal-baseline.csv", reference_date)))

write.csv(full, file = paste0(hub_outputpath,sprintf("%s-FluSight-seasonal-baseline.csv", reference_date)))

