# Flusight-baseline

The code in this repository creates the `Flusight-baseline` model and `Flusight-seasonal-baseline` model in the [Flusight-forecast-data](https://github.com/cdcepi/Flusight-forecast-data) repository.

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

The script `flusight-baseline.R` creates the Flusight-baseline forecasts. It can be run, for example, via

```
Rscript flusight-baseline.R
```

The script `flusight-seasonal-baseline.R` creates the Flusight-seasonal-baseline forecasts.

Baseline forecasts will be saved in the directory `weekly-submission/forecasts/Flusight-baseline`, and forecast plots will be saved in `weekly-submission/plots/Flusight-baseline`.

Seasonal baseline forecasts will be saved in the directory `weekly-submission/forecasts/Flusight-seasonal-baseline`.