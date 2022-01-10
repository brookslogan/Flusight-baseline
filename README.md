# Flusight-baseline

The code in this repository creates the `Flusight-baseline` model in the [Flusight-forecast-data](https://github.com/cdcepi/Flusight-forecast-data) repository.

## Dependencies

The baseline model is implemented in R, and has the following package dependencies:

 - Packages available on CRAN:
    - `purrr`
    - `dplyr`
    - `readr`
    - `tidyr`
    - `lubridate`
    - `ggplot2`
    - `ggforce`
 - Other packages, to be installed via GitHub:
    - `covidHubUtils`: See this package's [documentation page](http://reichlab.io/covidHubUtils/) for installation instructions.
    - `simplets`: Install from GitHub using the command `devtools::install_github("reichlab/simplets")`

## Running the baseline model

The script `baseline.R` creates the baseline forecasts. It can be run, for example, via

```
Rscript baseline.R
```

The script `fit_baseline_one_location.R` contains some helper functions and does not need to be run directly.

Baseline forecasts will be saved in the directory `weekly-submission/forecasts/Flusight-baseline`, and forecast plots will be saved in `weekly-submission/plots/Flusight-baseline`.
