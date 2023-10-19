# # The cdc_baseline_forecaster functionality is still in the v0.0.6 development
# # branch at the time of writing. Before it's merged, we'll work off of a
# # particular commit in the epipredict repo:
# pak::pkg_install("cmu-delphi/epipredict@3b809894d00c52ec9e166f26dc6f1c55c671b601")
#remotes::install_github("cmu-delphi/epipredict", ref="3b809894d00c52ec9e166f26dc6f1c55c671b601", dependencies = T)

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(checkmate)
library(cli)
library(epidatr)
library(epiprocess)
library(epipredict)
library(ggplot2)
library(plotly)

##############################
## Configuration parameters ##
##############################
output_dirpath <- "C:/Users/rpe5/Desktop/GitHub/Flusight-baseline/weekly-submission/forecasts/Flusight-baseline/"

######################
## Helper functions ##
######################

#' Return `date` if it has the desired weekday, else the next date that does
#' @param date `Date` vector
#' @param ltwday integerish vector; of weekday code(s), following POSIXlt
#'   encoding (not `lubridate`'s), but allowing either 0 or 7 to represent
#'   Sunday.
#' @return `Date` object
curr_else_next_date_with_ltwday <- function(date, ltwday) {
  assert_class(date, "Date")
  assert_integerish(ltwday, lower = 0L, upper = 7L)
  #
  date + (ltwday - as.POSIXlt(date)$wday) %% 7L
}

location_to_abbr <- function(location) {
  dictionary <-
    state_census %>%
    mutate(fips = sprintf("%02d", fips)) %>%
    transmute(
      location = case_match(fips, "00" ~ "US", .default = fips),
      abbr
    )
  dictionary$abbr[match(location, dictionary$location)]
}

###############################
## Fetch, prepare input data ##
###############################

# Get target data from cdcepi/FluSight-forecast-hub@main:
target_tbl <- readr::read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv",
  col_types = cols_only(
    date = col_date(format = ""),
    location = col_character(),
    location_name = col_character(),
    value = col_double(),
    weekly_rate = col_double()
  )
)

target_edf <- target_tbl %>%
  transmute(
    geo_value = location_to_abbr(location),
    time_value = .data$date,
    weekly_count = .data$value
  ) %>%
  as_epi_df()

# Implied date settings:
forecast_as_of_date <- Sys.Date()
reference_date <- curr_else_next_date_with_ltwday(forecast_as_of_date, 6L) # Saturday

# Validation:
desired_max_time_value <- reference_date - 14L

excess_latency_tbl <- target_edf %>%
  drop_na(weekly_count) %>%
  group_by(geo_value) %>%
  summarize(
    max_time_value = max(time_value),
    .groups = "drop"
  ) %>%
  mutate(
    excess_latency =
      pmax(
        as.integer(desired_max_time_value - max_time_value) %/% 7L,
        0L
      ),
    has_excess_latency = excess_latency > 0L
  )
excess_latency_small_tbl <- excess_latency_tbl %>%
  filter(has_excess_latency)

prop_locs_overlatent_err_thresh <- 0.20
prop_locs_overlatent <- mean(excess_latency_tbl$has_excess_latency)
if (prop_locs_overlatent > prop_locs_overlatent_err_thresh) {
  cli_abort("
    More than {100*prop_locs_overlatent_err_thresh}% of locations have excess
    latency. The reference date is {reference_date} so we desire observations at
    least through {desired_max_time_value}. However,
    {nrow(excess_latency_small_tbl)} location{?s} had excess latency and did not
    have reporting through that date: {excess_latency_small_tbl$geo_value}.
  ")
} else if (prop_locs_overlatent > 0) {
  cli_abort("
    Some locations have excess latency. The reference date is {reference_date}
    so we desire observations at least through {desired_max_time_value}.
    However, {nrow(excess_latency_small_tbl)} location{?s} had excess latency
    and did not have reporting through that date:
    {excess_latency_small_tbl$geo_value}.
  ")
}

######################
## Prepare baseline ##
######################

# For reproducibility, run with a particular RNG configuration. Make seed the
# same for all runs for the same `reference_date`, but different for different
# `reference_date`s. (It's probably not necessary to change seeds between
# `reference_date`s, though, since we use a large number of simulations so even
# if we sample the same quantile level trajectories, it won't be noticeable. The
# `%% 1e9` is also not necessary unless more seed-setting dependencies are added
# that would take us beyond R's integer max value.)
rng_seed <- as.integer((59460707 + as.numeric(reference_date)) %% 2e9)
withr::with_rng_version("4.0.0", withr::with_seed(rng_seed, {
  fcst <- cdc_baseline_forecaster(
    target_edf %>%
      # Match the start time_value filtering used in the 2022-23 baseline:
      filter(time_value >= as.Date("2021-12-04")) %>%
      # Don't use interim/preliminary data past the `desired_max_time_value`:
      filter(time_value <= desired_max_time_value),
    "weekly_count",
    cdc_baseline_args_list(
      # The `aheads` are specified relative to the most recent available
      # `time_value` available. Since our `data_frequency` is 1 week (the
      # default), the aheads are in terms of weeks.
      aheads = 1:5,
      nsims = 1e5,
      # (defaults for everything else)
    )
  )
  preds <- fcst$predictions %>%
    # epipredict infers a "`forecast_date`" equal to, and indexes aheads
    # relative to, the max `time_value` available, which is off from the
    # labeling we want due to data latency, but gives us the desired model and
    # `target_dates`. Instead, let the `forecast_date` be the `reference_date`
    # and index aheads relative to it:
    mutate(
      forecast_date = .env$reference_date,
      ahead = as.integer(.data$target_date - .env$reference_date) %/% 7L
    )
}))

##########
## Plot ##
##########

preds_wide <- pivot_quantiles_wider(preds, .pred_distn)
plot_states <- sort(unique(target_edf$geo_value))
plot_ncol <- 3L
plt <-
  preds_wide %>%
  filter(geo_value %in% plot_states) %>%
  mutate(geo_value = factor(geo_value, plot_states)) %>%
  arrange(geo_value) %>%
  ggplot(aes(target_date)) +
  geom_ribbon(aes(ymin = `0.1`, ymax = `0.9`), fill = blues9[3]) +
  geom_ribbon(aes(ymin = `0.25`, ymax = `0.75`), fill = blues9[6]) +
  geom_line(aes(y = .pred), color = "orange") +
  geom_line(
    data = target_edf %>%
      filter(geo_value %in% plot_states) %>%
      mutate(geo_value = factor(geo_value, plot_states)) %>%
      arrange(geo_value),
    aes(x = time_value, y = weekly_count)
  ) +
  scale_x_date(limits = c(reference_date - 120, reference_date + 30)) +
  labs(x = "Date", y = "Weekly admissions") +
  facet_wrap(~geo_value, scales = "free_y", ncol = plot_ncol) +
  # (as.numeric is workaround for
  # https://stackoverflow.com/questions/55150087/ggplotly-fails-with-geom-vline-with-xintercept-date-value)
  geom_vline(xintercept = as.numeric(desired_max_time_value), linetype = "dotted") +
  geom_vline(xintercept = as.numeric(forecast_as_of_date), linetype = "dotdash") +
  geom_vline(xintercept = as.numeric(reference_date), linetype = "dashed") +
  theme_bw()
ggplotly(plt, height = 400 * length(plot_states) / plot_ncol)

###################
## Format, write ##
###################

preds_formatted <- preds %>%
  flusight_hub_formatter(
    target = "wk inc flu hosp",
    output_type = "quantile"
  ) %>%
  drop_na(output_type_id) %>%
  arrange(target, horizon, location) %>%
  select(
    reference_date, horizon, target, target_end_date, location,
    output_type, output_type_id, value
  )

if (!dir.exists(output_dirpath)) {
  dir.create(output_dirpath, recursive = TRUE)
}
preds_formatted %>%
  write_csv(file.path(
    output_dirpath,
    sprintf("%s-FluSight-baseline.csv", reference_date)
  ))

