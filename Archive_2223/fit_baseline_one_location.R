#' Fit a baseline model  for one location
#'
#' Get quantile function
#'
#' @param  predictions baseline predictions
#' @param  taus probability levels
#'
get_quantiles_df <- function(predictions, taus) {
  n <- length(taus)
  purrr::map_dfr(
    1:4,
    function(h) {
      data.frame(
        h = rep(h, n),
        quantile = taus,
        value = ifelse(taus < 0.5,
          floor(quantile(predictions[, h], probs = taus)),
          ceiling(quantile(predictions[, h], probs = taus)))
      )
    })
}

#'
#' Get predictions
#'
#' @param  location_data data frame containing flu hospitalizations for a single location
#'   after outlier correction.
#' @param  response_var a value column after outlier detection and correction.
#' @param  transformation can be either "none" or "sqrt" or both.
#' @param  symmetrize can be either `TRUE` or `FALSE` or both.
#' @param  window_size a value or a vector of values of window size.
#' @param  h_adjust daily horizon adjustment for aggregation
#'
#' @return data frame of a baseline forecast for one location
get_baseline_predictions <- function(location_data,
                                     response_var,
                                     transformation,
                                     symmetrize,
                                     window_size,
                                     h_adjust,
                                     taus) {
  # fit
  baseline_fit <- fit_simple_ts(
    y = location_data[[response_var]],
    ts_frequency = 1,
    model = 'quantile_baseline',
    transformation = transformation,
    transform_offset = ifelse(transformation == "none", 0, 1),
    d = 0,
    D = 0,
    symmetrize = symmetrize,
    window_size = window_size
  )
  
  # predict
  weekly_predictions <-
    predict(baseline_fit, nsim = 100000, horizon = 4)
  
  # truncate to non-negative
  weekly_predictions <- pmax(weekly_predictions, 0)
  
  # extract predictive quantiles, intervals, and medians
  quantiles_df <- get_quantiles_df(weekly_predictions, taus)
  
  return(tibble(quantiles_df = list(quantiles_df)))
}

#' fit baseline to one location
#'
#' @param reference_date the date of the Saturday relative to which week-ahead targets are defined
#' @param location_data data frame containing flu hospitalizations for a single location. Must contain
#'   geo_value, time_value, and value columns.
#' @param transformation can be either "none" or "sqrt" or  both.
#' @param symmetrize can be either `TRUE` or `FALSE` or both.
#' @param window_size a value or a vector of values of window size.
#' @param taus probability levels
#'
#' @return data frame of a baseline forecast for one location
fit_baseline_one_location <- function(reference_date,
                                      location_data,
                                      transformation,
                                      symmetrize,
                                      window_size,
                                      taus) {
  predictions <- get_baseline_predictions(
    transformation = "none",
    symmetrize = TRUE,
    window_size = nrow(location_data),
    location_data = location_data,
    response_var = "value",
    taus = taus
  )

  # extract quantile forecasts
  quantiles_df <- predictions %>%
    tidyr::unnest(quantiles_df) %>%
    dplyr::transmute(
      forecast_date = as.character(forecast_date),
      target = paste0(h, " wk ahead inc flu hosp"),
      target_end_date = as.character(reference_date + 7L * h),
      location = unique(location_data$location),
      type = 'quantile',
      quantile = quantile,
      value = value,
      model = paste(
        "baseline",
        transformation,
        ifelse(symmetrize, "sym", "nonsym"),
        window_size,
        sep = "_"
      )
    )

  # add point estimates
  quantiles_df <- quantiles_df  %>%
    dplyr::bind_rows(
      .,
      quantiles_df %>%
        dplyr::filter(quantile == 0.5) %>%
        mutate(type = 'point',
               quantile = NA_real_)
    )

  return(quantiles_df)
}
