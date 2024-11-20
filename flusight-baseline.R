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
library(lubridate)

##############################
## Configuration parameters ##
##############################
userid <- Sys.info()["user"]
output_dirpath <- paste0("C:/Users/",userid,"/Desktop/GitHub/Flusight-baseline/weekly-submission/forecasts/Flusight-baseline/")
cat_ouput_dir <- paste0("C:/Users/",userid,"/Desktop/GitHub/FluSight-forecast-hub/model-output/FluSight-equal_cat/")

  
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

location_abbr_dictionary <-
  state_census %>%
  # convert fips to character if using old version of `state_census`:
  mutate(fips = if (is.numeric(fips)) sprintf("%02d", fips) else fips) %>%
  transmute(
    location = case_match(fips, "00" ~ "US", .default = fips),
    abbr
  )

location_to_abbr <- function(location) {
  location_abbr_dictionary$abbr[match(location, location_abbr_dictionary$location)]
}

abbr_to_location <- function(abbr) {
  location_abbr_dictionary$location[match(abbr, location_abbr_dictionary$abbr)]
}

###############################
## Fetch, prepare input data ##
###############################

target_tbl_col_spec <- cols_only(
  date = col_date(format = ""),
  location = col_character(),
  location_name = col_character(),
  value = col_double(),
  weekly_rate = col_double()
)
# (Reading in below tables with this col_spec may produce a message about
# renaming `` -> `...1` referencing the unnamed column containing row "names"
# (numbers), but cols_only will immediately drop it.)

# Final version of old-form (<= 2024) reporting:
target_tbl_old_form <-
  target_tbl <- readr::read_csv(
    "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/04e884dce942dd3b8766aee3d8ff1c333b4fb6fa/target-data/target-hospital-admissions.csv",
    col_types = target_tbl_col_spec
  )

# Latest version of new-form (>=2024) reporting mirrored at cdcepi/FluSight-forecast-hub@main:
target_tbl_new_form <- readr::read_csv(
  "https://raw.githubusercontent.com/cdcepi/FluSight-forecast-hub/main/target-data/target-hospital-admissions.csv",
  col_types = target_tbl_col_spec
)

if (min(target_tbl_new_form$date) <= as.Date("2022-12-31")) {
  # The new target table includes old time values spanning back pretty far; we
  # don't need to fill in with old target table values.
  target_tbl <- target_tbl_new_form
} else {
  # The new target table is missing a substantial time range from old-form
  # reporting. Fill in with old-form reporting where possible. Leave a time gap
  # between old-form and new-form reporting to prevent any jumps between the two
  # in the training set for the one-ahead model for the baseline.
  target_tbl <- bind_rows(
    target_tbl_old_form %>% filter(date <= min(target_tbl_new_form$time_value) - 14L),
    target_tbl_new_form
  )
}
# We'll also filter out some early time values below when training the model.

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
desired_max_time_value <- reference_date - 7L

# * that we're not running too late:
max_time_value <- max(target_edf$time_value)
if (max_time_value > desired_max_time_value) {
  cli_abort("
    The target data run through a max time value of {max_time_value},
    but we were expecting them to run only through {desired_max_time_value}
    in order to make predictions at forecast date {forecast_as_of_date},
    reference date {reference_date}.
  ")
}

# * that data's not running too late / we're not running too early:
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

imposed_min_time_value <- as.Date("2022-08-06") # 2022EW31 Sat
#
# ^ For seasons through 2023/2024, this was instead 2021-12-04. For 2024/2025,
# it has been updated to exclude the low activity during 2021/2022. EW31 was
# selected as a boundary between 2021/2022 and 2022/2023 to nearly-evenly divide
# up off-season weeks and to include the full 2022/2023 season ramp-up, though
# this also includes more flat off-season weeks.

pause_min_time_value <- as.Date("2024-04-27") + 7L # Sat
pause_max_time_value <- as.Date("2024-11-09") - 7L # Sat

n_output_trajectories <- 100L

# For reproducibility, run with a particular RNG configuration. Make seed the
# same for all runs for the same `reference_date`, but different for different
# `reference_date`s. (It's probably not necessary to change seeds between
# `reference_date`s, though, since we use a large number of simulations so even
# if we sample the same quantile level trajectories, it won't be noticeable. The
# `%% 1e9` is also not necessary unless more seed-setting dependencies are added
# that would take us beyond R's integer max value.)
rng_seed <- as.integer((59460707 + as.numeric(reference_date)) %% 2e9)
withr::with_rng_version("4.0.0", withr::with_seed(rng_seed, {
  # Temporary approach to grab trajectory samples from epipredict without
  # requiring epipredict update:
  subsamples_by_geo <- list()
  trace(epipredict:::propagate_samples, exit = quote({
    n <- 1L
    e <- rlang::caller_env(n)
    while(! ".data" %in% names(e)) {
      n <- n + 1L
      e <- rlang::caller_env(n)
      if (identical(e, globalenv()) || identical(e, emptyenv())) {
        cli::cli_abort("Failed to find the target geo_value to attach to this trajectory sample.")
      }
    }
    target_geo <- e$.data$geo_value
    sample_by_horizon <- res
    # sample_by_horizon is sorted by horizon 0 draws, so we can't just take the
    # first n_output_trajectories of them; subsample instead:
    selected_trajectory_inds <- sample.int(length(sample_by_horizon[[1L]]), n_output_trajectories)
    subsample_by_horizon <- lapply(sample_by_horizon, `[`, selected_trajectory_inds)
    # Ensure non-negative:
    subsample_by_horizon <- lapply(subsample_by_horizon, pmax, 0L)
    # Prepare sample ids; the ith draw for each of the horizons belong to the
    # same sample, so they should have the same sample id; sampling is performed
    # separately for each geo, so different geos should have different sets of
    # ids:
    subsample_ids_every_horizon <- paste0(target_geo, "_s", seq_len(n_output_trajectories))
    subsample <- tibble(
      geo_value = target_geo,
      horizon = seq_along(subsample_by_horizon) - 1L,
      output_type_id = rep(list(subsample_ids_every_horizon), length(horizon)),
      value = subsample_by_horizon
    ) %>%
      unchop(c(output_type_id, value))
    .GlobalEnv[["subsamples_by_geo"]] <- c(.GlobalEnv[["subsamples_by_geo"]], list(subsample))
  }))
  #
  # Forecasts for all but the -1 horizon, in `epipredict`'s forecast output
  # format. We will want to edit some of the labeling and add horizon -1, so we
  # won't use this directly.
  fcst <- cdc_baseline_forecaster(
    target_edf %>%
      filter(time_value >= imposed_min_time_value) %>%
      filter(!between(time_value, pause_min_time_value, pause_max_time_value)) %>%
      # Don't use interim/preliminary data past the `desired_max_time_value`
      # (shouldn't do anything if we raised an error earlier on about
      # unexpectedly low latency):
      filter(time_value <= desired_max_time_value),
    "weekly_count",
    cdc_baseline_args_list(
      # The `aheads` are specified relative to the most recent available
      # `time_value` available. Since our `data_frequency` is 1 week (the
      # default), the aheads are in terms of weeks.
      aheads = 1:4,
      nsims = 1e5,
      # (defaults for everything else)
    )
  )
  # Extract the predictions in `epipredict` format, and add horizon -1
  # predictions:
  preds <- fcst$predictions %>%
    # epipredict infers a "`forecast_date`" equal to, and indexes aheads
    # relative to, the max `time_value` available, which is off from the
    # labeling we want due to data latency, but gives us the desired model and
    # `target_dates`. Instead, let the `forecast_date` be the `reference_date`
    # and index aheads relative to it:
    mutate(
      forecast_date = .env$reference_date,
      ahead = as.integer(.data$target_date - .env$reference_date) %/% 7L
    ) %>%
    bind_rows(
      # Prepare -1 horizon predictions:
      target_edf %>%
        # Pretend that excess latency, either in the form of missing rows or
        # NAs, doesn't exist; the last available week will be treated as if it
        # ended on `desired_max_time_value`:
        drop_na(weekly_count) %>%
        slice_max(time_value) %>%
        transmute(
          # Like in the preceding rows of `preds`, we will let `forecast_date`
          # be the `reference_date` and index aheads relative to it:
          forecast_date = .env$reference_date,
          target_date = .env$reference_date - 7L,
          ahead = -1L,
          geo_value,
          .pred = weekly_count,
          # Degenerate (deterministic) distributions:
          .pred_distn = dist_quantiles(
            values = map(
              weekly_count, rep,
              length(cdc_baseline_args_list()$quantile_levels)
            ),
            quantile_levels = cdc_baseline_args_list()$quantile_levels
          )
        )
    )
}))

subsamples <- subsamples_by_geo %>%
  bind_rows() %>%
  mutate(reference_date = .env$reference_date,
         target_date = reference_date + 7L * horizon) %>%
  mutate(value = round(value)) # currently inconsistent with quantile preds not rounding

##########
## Plot ##
##########

preds_wide <- pivot_quantiles_wider(preds, .pred_distn)
plot_states <- sort(unique(target_edf$geo_value))
plot_ncol <- 3L
plot_ntraj <- 10L
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
  geom_line(
    data = subsamples %>% slice_head(n = plot_ntraj, by = c(geo_value, horizon)),
    aes(y = value, group = output_type_id),
    alpha = 0.1
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

quantile_preds_formatted <- preds %>%
  flusight_hub_formatter(
    target = "wk inc flu hosp",
    output_type = "quantile"
  ) %>%
  drop_na(output_type_id) %>%
  arrange(target, horizon, location) %>%
  dplyr::mutate(
    value = ifelse(output_type_id < 0.5, floor(value), ceiling(value))  # Round value based on output_type_id
  ) %>%
  dplyr::mutate(
    output_type_id = as.character(output_type_id)
  ) %>%
  dplyr::select(
    reference_date, horizon, target, target_end_date, location,
    output_type, output_type_id, value
  )

sample_preds_formatted <- subsamples %>%
  transmute(
    reference_date,
    horizon,
    target = "wk inc flu hosp",
    target_end_date = target_date,
    location = abbr_to_location(geo_value),
    output_type = "sample",
    output_type_id,
    value
  )

preds_formatted <- bind_rows(list(
  quantile_preds_formatted,
  sample_preds_formatted
))

if (!dir.exists(output_dirpath)) {
  dir.create(output_dirpath, recursive = TRUE)
}
preds_formatted %>%
  write_csv(file.path(
    output_dirpath,
    sprintf("%s-FluSight-baseline.csv", reference_date)
  ))



###########################################################
## Prepare flat baseline for categorical trend forecasts ##
###########################################################

# This code uses a template data file for with equal probabilities for each of the five categories
# and then updates the dates to correspond with a submission with the upcoming Saturday as a reference
# date. 

flat_cat_template <- read.csv(paste0("C:/Users/",userid,"/Desktop/GitHub/FluSight-forecast-hub/model-output/FluSight-equal_cat/2024-01-06-FluSight-equal_cat.csv"),header=T)

# Find the next Saturday for updating reference_date
(nsat<-today()+days(7-wday(today())))
# Update reference_date and target_end_dates for horizon = 0,1,2,3 
next_dates<-flat_cat_template %>% 
  mutate(reference_date=nsat,
         target_end_date=reference_date+weeks(horizon))
write.csv(next_dates,file=paste0(cat_ouput_dir,"/",nsat,"-FluSight-equal_cat.csv"),row.names=FALSE)

## Backup for manually entering a Saturday (e.g., past reference_date needed)
# saturday<-ymd("2023-10-21")
# new_date<-template %>% 
#   mutate(reference_date=saturday,
#          target_end_date=reference_date+weeks(template$horizon))
# write.csv(new_date,file=paste0(outputdir,"/",saturday,"-FluSight-equal_cat.csv"),row.names=FALSE)
