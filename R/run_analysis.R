## packages
pacman::p_load(
  meta, metafor, lmtest, scimple, tidyverse, rio, broom, here,
  glue, magrittr, warp, binom, nnet, effects, MultinomialCI,
  RColorBrewer, patchwork, ggtext, purrr, ggrepel, splitstackshape,
  rfextras, DescTools
)

## functions
source(here("R/functions.R"))

## get gisaid genomic epi data
gis <- get_gis()

## ## get epinow
## epinow <- get_epinow()

## this will run the analysis with the 4 most common lineages
results <- get_results(
  ## provide data
  gis = gis,
  ## do you want to re-run the analysis?
  update = TRUE,
  ## what is the date of GISAID database?
  db_date = Sys.Date(),
  ## how many variants do you want to investigate per country?
  n_var = 4,
  ## when calculating which variants to look at, we take
  ## the most common from the last x days - how many days
  ## should we look back?
  last_n_days = 50,
  ## if you want to look at specific variants instead of
  ## the most common recent ones, specify them here
  ## variants = c("B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "B.1.617.1",
  ##              "P.2", "P.3", "B.1.525", "B.1.526", "B.1.427/B.1.429"),
  ## minimum number of variant sequences to include for analysis
  min_n_var = 15,
  ## how many days to forecast for
  days_forecast = 14
)

## this will run the analysis with only VOI/VOC
results <- get_results(
  ## provide data
  gis = gis,
  ## do you want to re-run the analysis?
  update = TRUE,
  ## what is the date of GISAID database?
  db_date = Sys.Date(),
  ## how many variants do you want to investigate per country?
  n_var = 4,
  ## when calculating which variants to look at, we take
  ## the most common from the last x days - how many days
  ## should we look back?2
  last_n_days = 50,
  ## if you want to look at specific variants instead of
  ## the most common recent ones, specify them here
  variants = c("B.1.1.7", "B.1.351", "P.1", "B.1.617.2", "B.1.617.1",
               "P.2", "P.3", "B.1.525", "B.1.526", "B.1.427/B.1.429"),
  ## minimum number of variant sequences to include for analysis
  min_n_var = 15,
  ## how many days to forecast for
  days_forecast = 14
)

## look at combined trajectories
results$plot_combined[[3]]

## look at split trajectories
results$plot_split[[3]]

## look at split trajectories
results$plot_rt[[29]]

## look at raw data
vis_raw(
  ## provide data
  gis = gis,
  ## which countries to look at?
  countries = c("India", "United Kingdom", "Brazil"),
  ## how many weeks to group data by
  weeks = 1,
  ## should errorbars be shown?
  errorbars = TRUE,
  ## how many variants should be looked at?
  n_variants = 5,
  ## when calculating which variants to look at, we take
  ## the most common from the last x days - how many days
  ## should we look back?
  last_n_days = 50,
  ## or manually specify variants
  variants = c("B.1.617.1", "B.1.617.2", "B.1.1.7", "P.1"),
  ## what range of dates should we look at?
  date_range = c(as.Date("2020-06-01"), Sys.Date()),
  ## specify brewer palette
  palette = "Set1",
  ## set font size
  base_size = 16
)
