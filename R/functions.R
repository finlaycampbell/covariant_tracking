## pull covariate data
pull_cov_n501y <- function() {

  dat <- jsonlite::fromJSON(
    url(paste0("https://raw.githubusercontent.com/hodcroftlab/",
               "covariants/master/cluster_tables/S.N501_data.json"))
  )

  parse_subs <- function(data, name) {
    data <- as_tibble(data)
    data <- mutate(data, country = name) %>%
      relocate(country, .before = 1)
    return(data)
  }

  covariants_data <- purrr::map2_df(dat, names(dat), ~parse_subs(.x, .y)) %>%
    mutate(week = lubridate::ymd(week),
           prop_n501 = cluster_sequences / total_sequences) %>%
    mutate(
      country = case_when(
        country == "USA" ~ "United States of America",
        TRUE ~ country
      )
    )

  right_join(
    phifunc::pull_pop_data() %>%
    select(report_country, iso3, who_region),
    covariants_data,
    by = c("report_country" = "country")
  )

}

## fit by growth rate
fit_growth <- function(df) {
  nls(
    y ~ start/(start + (1 - start)*exp(-t*diff)),
    data = df,
    start = c(start = df$y[1], diff = 0.1)
  )
}

## fit by reproduction number
fit_r <- function(df, gen_time = 5.2) {
  nls(
    y ~ start/(start + (1 - start)*ratio^(-t/gen_time)),
    data = df,
    start = c(start = df$y[1], ratio = 1)
  )
}

## forecasting growth model
forecast_growth <- function(tidied, data,
                            t_start = as.Date("2020-10-01"),
                            t_end = Sys.Date() + 45) {

  root <- min(data$week)
  dates <- seq.Date(t_start, t_end, by = 1)
  ind <- as.numeric(dates - root)
  ## project across upper mid and lower

  pmap_dfc(
    tidied,
    ~ tibble(
      var = ..2/(..2 + (1 - ..2)*exp(-(ind)*..3))
    )
  ) %>%
    setNames(tidied$conf) %>%
    mutate(date = dates)

}

## forecasting r model
forecast_r <- function(tidied, data,
                       gen_time = 5.2,
                       t_start = as.Date("2020-10-01"),
                       t_end = Sys.Date() + 45) {


  root <- min(data$week)
  dates <- seq.Date(t_start, t_end, by = 1)
  ind <- as.numeric(dates - root)
  ## project across upper mid and lower

  pmap_dfc(
    tidied,
    ~ tibble(
      var = ..2/(..2 + (1 - ..2)*..3^(-ind/gen_time))
    )
  ) %>%
    setNames(tidied$conf) %>%
    mutate(date = dates)

}

## forecasting r model
forecast_glm_manual <- function(tidied, data,
                                t_start = as.Date("2020-10-01"),
                                t_end = Sys.Date() + 45) {
  root <- min(data$date)
  dates <- seq.Date(t_start, t_end, by = 1)
  ind <- as.numeric(dates - root)
  ## project across upper mid and lower
  imap_dfr(
    select(tidied, -term),
    ~ tibble(
      val =  .x[1]/(.x[1] + (1 - .x[1])*exp(-.x[2]*ind)),
      var = .y,
      date = dates
    )
  ) %>%
    pivot_wider(names_from = var, values_from = val)
}

## weirdly necessary global variable mapping for as.lm function
mymap <- function(x, f, ...) {
  x <<- x
  out <- lst()
  for(i in seq_along(x)) out[[i]] <- f(x[[i]], ...)
  return(out)
}

## plot saving function
save_plot <- function(p, file,
                      width = 20.16, height = 13.35,
                      units = 'cm',
                      dpi = 600,
                      folder = "figures",
                      ...) {

  if(is.null(p)) return(NULL)

  if(!file.exists(here(folder))) dir.create(here(folder))

  ggsave(
    filename = here(folder, file),
    plot = p,
    width = width,
    height = height,
    units = units,
    dpi = dpi,
    ...
    )

}

## visualise fits
vis_nls <- function(fits, min_rsq = 0.9) {

  fits %<>%
    filter(
      growth_rsq > min_rsq,
      map_dbl(growth_tidied, pluck, "lower", 2, .default = NA) > 0
    ) %>%
    group_by(report_country) %>%
    mutate(
      est = map_chr(
        growth_tidied,
        possibly(
          ~ glue("{sprintf('%.2f', .x$mid[2])} [{sprintf('%.2f', .x$lower[2])}",
                        " - {sprintf('%.2f', .x$upper[2])}]"), NA)
      ),
      lab = glue("{report_country}\n{est[which.max(growth_rsq)]}")
    )

  df1 <- unnest(fits, data)
  df2 <- unnest(fits, growth_forecast) %>%
    mutate(
      across(c(lower, upper, mid), ~ replace(.x, .x < 0, 0)),
      across(c(lower, upper, mid), ~ replace(.x, .x > 1, 1))
    )

  ggplot() +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    geom_line(
      data = df2,
      aes(date, y = mid)
    ) +
    geom_point(
      data = df1,
      aes(week, y = prop_n501),
      size = 0.8,
      alpha = 0.5
    ) +
    labs(x = NULL, y = "Proportion B.1.1.7") +
    facet_wrap(~ lab) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

}

## get palette for strains
get_var_pal <- function(labs, pal = "Dark2", pkg = c("brewer", "tableau"), other_col = NULL) {

  pkg <- match.arg(pkg)

  labs <- if(is.null(other_col)) head(labs, -1) else head(labs, -1)[-1]

  cols <- if(pkg == "brewer") get_brewer(length(labs), pal) else get_tableau(length(labs), pal)

  names(cols) <- names(labs)

  if(!is.null(other_col)) cols <- c(str1 = other_col, cols)

  cols %<>% modify_if(~ .x == "#BAB0AC", ~ "black") %>%
    modify_if(~ .x == "#FF9DA7", ~ "#BAB0AC")

  return(cols)

}

## relabel variant definitions
get_var_lab <- function(vars, md = FALSE, swap = TRUE, remove = FALSE) {

  mtch <- c(
    B.1.1.7 = "Alpha",
    B.1.351 = "Beta",
    P.1 = "Gamma",
    B.1.617.2 = "Delta",
    B.1.617.1 = "Kappa",
    P.2 = "Zeta",
    P.3 = "Theta",
    B.1.525 = "Eta",
    B.1.526 = "Iota",
    `B.1.427/B.1.429` = "Epsilon"
  )


  if(!remove)
    vars %<>% modify_if(
      ~ .x %in% names(mtch),
      if(swap) {if(md) ~ glue("{.x}<br>({mtch[.x]})") else ~ glue("{.x}\n({mtch[.x]})")}
      else {if(md) ~ glue("{mtch[.x]}<br>({.x})") else ~ glue("{mtch[.x]}\n({.x})")}
    )
  else
    vars %<>% modify_if(
      ~ .x %in% names(mtch),
      ~ mtch[.x]
    )

  return(vars)

}

## get color brewer colors
get_brewer <- function(n, pal = "Dark2") {

  nmax <- brewer.pal.info$maxcolors[which(rownames(brewer.pal.info) == pal)]
  if(n <= nmax) {
    cols <- brewer.pal(nmax, pal)[seq_len(n)]
  } else {
    cols <- colorRampPalette(brewer.pal(nmax, pal))(n)
  }

  return(cols)

}

## get tableau colors
get_tableau <- function(n, pal = "Tableau 10") {

  raw <- ggthemes::ggthemes_data[["tableau"]][["color-palettes"]][[c("regular")]][[pal]]

  if(pal == "Tableau 10") raw$value <- raw$value[c(1:4, 6, 5, 7:10)]

  if(n <= nrow(raw)) {
    cols <- raw$value[seq_len(n)]
  } else {
    cols <- colorRampPalette(raw$value)(n)
  }

  return(cols)

}

## visualise fits
vis_mglm <- function(fits,
                     countries = unique(fits$report_country),
                     weeks = 2,
                     log = FALSE,
                     date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                     top_n = 30,
                     threshold = 0.001) {

  ## stack values across strainsg
  stack_values <- function(df, n) {
    for(i in 2:n) {
      vals <- df[[paste0("mid_str", i - 1)]]
      df %<>% mutate(
        across(
          paste0(c("lower_str", "mid_str", "upper_str"), as.character(i)),
          ~ .x + vals
        )
      )
    }
    return(df)
  }

  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>%
    filter(!map_lgl(forecast, is.null) & report_country %in% countries) ## %>%
    ## mutate(report_country = paste0(report_country, " - ", round(pseudo, 2)))

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      slice(1:top_n)
  }

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %T>%
    {n_var <<- length(levels(.$y))} %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  ## remove these points in the future
  to_na <- which(df1 < threshold, arr.ind = TRUE)

  ## stack values
  df1 %<>% stack_values(n_var)

  ## ## remove small values
  if(nrow(to_na) > 0) {
    for(i in seq_len(nrow(to_na))) df1[to_na[i,1], to_na[i,2]] <- NA
  }

  ## transform
  df1 %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## forecast results
  df2 <- unnest(fits, forecast) %>%
    filter(date > date_range[1] & date < date_range[2])

  ## stack forecast results and transform
  df2_adj <- stack_values(df2, n_var) %>% piv_long() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))
  df2 %<>% piv_long() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## adjust for log results
  if(log) {
    df2 %<>% filter(mid != 0)
    df1 %<>% filter(mid != 0)
  }

  df1$strain %<>% droplevels()
  df2 %<>%
    filter(as.character(strain) %in% levels(df1$strain)) %>%
    mutate(strain = droplevels(strain))
  df2_adj$strain %<>% droplevels()

  cols <- get_var_pal(fits$labs[[1]])

  ggplot() +
    ## geom_vline(xintercept = Sys.Date(), linetype = 2) +
    geom_area(
      data = df2,
      aes(date, y = mid, fill = fct_rev(strain)),
      alpha = 0.4,
      color = NA
    ) +
    geom_ribbon(
      data = df2_adj,
      aes(date, ymin = lower, ymax = upper, group = strain),
      alpha = 1,
      color = "white",
      fill = "white"
    ) +
    geom_line(
      data = df1,
      aes(week_one, y = mid, group = fct_rev(strain)),
      color = 'black',
      alpha = 0.5,
      size = 0.3
    ) +
    geom_point(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.7,
      alpha = 1
    ) +
    geom_errorbar(
      data = df1,
      aes(week_one, ymin = lower, ymax = upper, color = fct_rev(strain)),
      width = 0,
      size = 0.3
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0, 0),
      ## date_breaks = "2 months",
      limits = range(df2$date),
      date_labels = "%b\n%Y"
    ) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_fill_manual(
      values = cols,
      name = "Variant",
      labels = fits$labs[[1]],
    ) +
    scale_color_manual(
      values = cols,
      name = "Variant",
      labels = fits$labs[[1]],
    ) +
    guides(
      fill = guide_legend(reverse = TRUE),
      color = guide_legend(reverse = TRUE)
    ) +
    labs(x = NULL, y = "Proportion of sequences") +
    facet_wrap(~ report_country, ncol = 5) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom'
    )

}

## visualise fits
vis_ribbons <- function(fits,
                     countries = unique(fits$report_country),
                     weeks = 2,
                     log = FALSE,
                     date_range = c(as.Date("2020-09-01"), Sys.Date()),
                     top_n = 30,
                     days_forecast = 14,
                     ncol = 5,
                     base_size = 9,
                     axis_text_size = 6,
                     scaler = 1,
                     rem_most_recent = c("South Africa", "India", "United Kingdom"),
                     threshold = 0.001) {

  ## stack values across strainsg
  stack_values <- function(df, n) {
    for(i in 2:n) {
      vals <- df[[paste0("mid_str", i - 1)]]
      df %<>% mutate(
        across(
          paste0(c("lower_str", "mid_str", "upper_str"), as.character(i)),
          ~ .x + vals
        )
      )
    }
    return(df)
  }

  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  ## save labs
  labs <- fits$labs[[1]]

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>%
    filter(!map_lgl(forecast, is.null) & report_country %in% countries) ## %>%
    ## mutate(report_country = paste0(report_country, " - ", round(pseudo, 2)))

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      slice(1:top_n)
  }

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %T>%
    {n_var <<- length(levels(.$y))} %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  df1_adj <- df1

  ## remove these points in the future
  to_na <- which(df1_adj < threshold, arr.ind = TRUE)

  ## stack values
  df1_adj %<>% stack_values(n_var)

  ## ## remove small values
  if(nrow(to_na) > 0) {
    for(i in seq_len(nrow(to_na))) df1_adj[to_na[i,1], to_na[i,2]] <- NA
  }

  ## transform
  df1_adj %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## transform
  df1 %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## forecast results
  df2 <- unnest(fits, forecast) ## %>%
    ## filter(date > date_range[1] & date < date_range[2])

  ## stack forecast results and transform
  df2_adj <- stack_values(df2, n_var) %>% piv_long() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))
  df2 %<>% piv_long() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## adjust for log results
  if(log) {
    df2 %<>% filter(mid != 0)
    df1 %<>% filter(mid != 0)
  }

  df1$strain %<>% droplevels()
  df2 %<>%
    filter(as.character(strain) %in% levels(df1$strain)) %>%
    mutate(strain = droplevels(strain))
  df2_adj$strain %<>% droplevels()

  cols <- get_var_pal(fits$labs[[1]], pal = "Tableau 10", pkg = "tableau")

  txt <- highlight_voc(fits$labs[[1]], remove = TRUE)

  min_date <- df1 %>%
    group_by(report_country) %>%
    summarise(
      min_date = min(week_one),
      max_date = max(week_one) + days_forecast
    )

  df2 %<>%
    left_join(min_date, "report_country") %>%
    filter(date >= min_date & date <= max_date)

  df2_adj %<>%
    left_join(min_date, "report_country") %>%
    filter(date >= min_date & date <= max_date)

  df1 %<>%
    group_by(report_country) %>%
    filter(!(week_one == max(week_one) & report_country %in% rem_most_recent))

  df1_adj %<>%
    group_by(report_country) %>%
    filter(!(week_one == max(week_one) & report_country %in% rem_most_recent))

  ggplot() +
    ## geom_vline(xintercept = Sys.Date(), linetype = 2) +
    geom_area(
      data = df2,
      aes(date, y = mid, fill = fct_rev(strain), color = NA),
      alpha = 0.5,
      color = NA
    ) +
    geom_ribbon(
      data = df2_adj,
      aes(date, ymin = lower, ymax = upper, group = strain),
      alpha = 1,
      color = "white",
      fill = "white"
    ) +
    geom_area(
      data = df1,
      aes(week_one, y = mid, ymin = lower, ymax = upper,
          group = fct_rev(strain), color = fct_rev(strain), fill = fct_rev(strain)),
      size = 0.3 * scaler,
      alpha = 0.5
    ) +
    geom_point(
      data = df1_adj,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.4 * scaler,
      alpha = 1
    ) +
    geom_errorbar(
      data = df1_adj,
      aes(week_one, ymin = lower, ymax = upper, color = fct_rev(strain)),
      width = 0,
      size = 0.3 * scaler
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "2 month",
      limits = range(df2$date),
      date_labels = "%b\n%Y"
    ) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_fill_manual(
      values = cols,
      name = "Variant",
      labels = txt,
    ) +
    scale_color_manual(
      values = cols,
      name = "Variant",
      labels = txt,
    ) +
    guides(
      fill = guide_legend(reverse = TRUE),
      color = guide_legend(reverse = TRUE)
    ) +
    labs(x = NULL, y = "Proportion of sequences") +
    facet_wrap(~ report_country, ncol = ncol) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.text = element_markdown(),
      axis.text.x = element_text(size = axis_text_size)
    )

}

## visualise fits
vis_modelled <- function(fits,
                         countries = unique(fits$report_country),
                         weeks = 2,
                         log = FALSE,
                         date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                         top_n = 30,
                         threshold = 0.001) {

  ## stack values across strainsg
  stack_values <- function(df, n) {
    for(i in 2:n) {
      vals <- df[[paste0("mid_str", i - 1)]]
      df %<>% mutate(
        across(
          paste0(c("lower_str", "mid_str", "upper_str"), as.character(i)),
          ~ .x + vals
        )
      )
    }
    return(df)
  }

  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  ## save labs
  labs <- fits$labs[[1]]

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>%
    filter(!map_lgl(forecast, is.null) & report_country %in% countries) ## %>%
    ## mutate(report_country = paste0(report_country, " - ", round(pseudo, 2)))

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      slice(1:top_n)
  }

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %T>%
    {n_var <<- length(levels(.$y))} %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  df1_adj <- df1

  ## remove these points in the future
  to_na <- which(df1_adj < threshold, arr.ind = TRUE)

  ## stack values
  df1_adj %<>% stack_values(n_var)

  ## ## remove small values
  if(nrow(to_na) > 0) {
    for(i in seq_len(nrow(to_na))) df1_adj[to_na[i,1], to_na[i,2]] <- NA
  }

  ## transform
  df1_adj %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## transform
  df1 %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## forecast results
  df2 <- unnest(fits, forecast) %>%
    filter(date > date_range[1] & date < date_range[2])

  ## stack forecast results and transform
  df2_adj <- stack_values(df2, n_var) %>% piv_long() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))
  df2 %<>% piv_long() %>%
    mutate(strain = factor(strain, head(names(fits$labs[[1]]), -1)))

  ## adjust for log results
  if(log) {
    df2 %<>% filter(mid != 0)
    df1 %<>% filter(mid != 0)
  }

  df1$strain %<>% droplevels()
  df2 %<>%
    filter(as.character(strain) %in% levels(df1$strain)) %>%
    mutate(strain = droplevels(strain))
  df2_adj$strain %<>% droplevels()

  cols <- get_var_pal(fits$labs[[1]])

  txt <- modify_if(
    fits$labs[[1]],
    ~ .x %in% c("P.1", "B.1.1.7", "B.1.351", "B.1.617.1", "B.1.617.2", "B.1.617.3"),
    ~ paste0("**", .x, "**")
  )

  min_date <- df1 %>%
    group_by(report_country) %>%
    summarise(min_date = min(week_one))

  df2 %<>%
    left_join(min_date, "report_country") %>%
    filter(date >= min_date)

  df2_adj %<>%
    left_join(min_date, "report_country") %>%
    filter(date >= min_date)

  ggplot() +
    geom_area(
      data = df2,
      aes(date, y = mid, fill = fct_rev(strain), color = NA),
      alpha = 1,
      color = "black"
    ) +
    ## geom_ribbon(
    ##   data = df2_adj,
    ##   aes(date, ymin = lower, ymax = upper, group = strain),
    ##   alpha = 1,
    ##   color = "white",
    ##   fill = "white"
    ## ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0, 0),
      date_breaks = "2 month",
      limits = range(df2$date),
      date_labels = "%b\n%Y"
    ) +
    coord_cartesian(ylim=c(0, 1)) +
    scale_fill_manual(
      values = cols,
      name = "Variant",
      labels = txt,
    ) +
    scale_color_manual(
      values = cols,
      name = "Variant",
      labels = txt,
    ) +
    guides(
      fill = guide_legend(reverse = TRUE),
      color = guide_legend(reverse = TRUE)
    ) +
    labs(x = NULL, y = "Proportion of sequences") +
    {if(length(unique(df2$report_country)) > 1) facet_wrap(~ report_country, ncol = 5) else NULL} +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom',
      legend.text = element_markdown()
    )

}

## visualise fits
vis_raw_old <- function(gis,
                    countries,
                    n_variants = 5,
                    last_n_days = 50,
                    variants = NULL,
                    weeks = 2,
                    date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                    palette = "Dark2",
                    base_size = 9) {

  ## nest data
  fits <- nest(gis, data = -report_country) %>%
    filter(report_country %in% countries)

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(as.data.frame(prop.table(table(.$pango_lineage)))) %>%
    rename(strain = Var1, mid = Freq) %>%
    ungroup()

  if(is.null(variants)) {
    variants <- get_most_common(
      gis, countries, n_var = n_variants,
      date_range = c(Sys.Date() - last_n_days, Sys.Date())
    ) %$%
      c(pango_lineage, "Other")
  } else variants <- c(variants, "Other")

  df1 %<>%
    mutate(
      strain = as.character(strain),
      strain = replace(strain, !strain %in% variants, "Other"),
      strain = factor(strain, variants)
    ) %>%
    group_by(report_country, week_one, strain) %>%
    summarise(mid = sum(mid)) %>%
    complete(report_country, week_one, strain, fill = list(mid = 0))

  get_pal <- function(n, pal = "Dark2") {
    nmax <- brewer.pal.info$maxcolors[which(rownames(brewer.pal.info) == pal)]
    if(n <= nmax)
      cols <- brewer.pal(nmax, pal)[seq_len(n)]
    else
      cols <- colorRampPalette(brewer.pal(nmax, pal))(n)
    return(cols)
  }

  ggplot() +
    geom_area(
      data = df1,
      aes(week_one, y = mid, fill = strain),
      color = 'black'
    ) +
    scale_y_continuous(
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y"
    ) +
    scale_fill_manual(values = get_pal(length(unique(df1$strain)), palette)) +
    labs(x = NULL, y = "Proportion of sequences", fill = "Variant") +
    {if(length(countries) > 1) facet_wrap(~ report_country)} +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom'
    )

}

## visualise fits
vis_raw <- function(gis,
                    countries,
                    n_variants = 5,
                    last_n_days = 50,
                    variants = NULL,
                    weeks = 2,
                    date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                    palette = "Dark2",
                    errorbars = TRUE,
                    threshold = 0.001,
                    base_size = 9) {


  ## stack values across strainsg
  stack_values <- function(df, n) {
    for(i in 2:n) {
      vals <- df[[paste0("mid_str", i - 1)]]
      df %<>% mutate(
        across(
          paste0(c("lower_str", "mid_str", "upper_str"), as.character(i)),
          ~ .x + vals
        )
      )
    }
    return(df)
  }

  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  if(is.null(variants)) {
    variants <- get_most_common(
      gis, countries, n_var = n_variants,
      date_range = c(Sys.Date() - last_n_days, Sys.Date())
    ) %$%
      c("Other", pango_lineage)
  } else variants <- c("Other", variants)

  ## nest data
  fits <- nest(gis, data = -report_country) %>%
    filter(report_country %in% countries)

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    mutate(
      dist_one = warp_distance(date, "week", every = weeks),
      y = factor(
        paste0("str", replace_na(match(pango_lineage, variants[-1]) + 1, 1)),
        levels = paste0("str", seq_len(length(variants)))
      )
    ) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  df1_adj <- df1

  ## remove these points in the future
  to_na <- which(df1_adj < threshold, arr.ind = TRUE)

  ## stack values
  df1_adj %<>% stack_values(length(variants))

  ## ## remove small values
  if(nrow(to_na) > 0) {
    for(i in seq_len(nrow(to_na))) df1_adj[to_na[i,1], to_na[i,2]] <- NA
  }

  ## stack values
  df1_adj %<>%
    piv_long() %>%
    drop_na() %>%
    mutate(
      strain = factor(strain, paste0("str", seq_len(length(variants)))),
      lower = replace(lower, lower < 0, 0),
      upper = replace(upper, upper > 1, 1)
    ) %>%
    group_by(report_country) %>%
    complete(report_country, week_one, strain, fill = list(lower = 0, mid = 0, upper = 0)) %>%
    ungroup()

  ## transform
  df1 %<>% piv_long() %>%
    drop_na() %>%
    mutate(strain = factor(strain, paste0("str", seq_len(length(variants))))) %>%
    select(-lower, -upper) %>%
    group_by(report_country) %>%
    complete(report_country, week_one, strain, fill = list(mid = 0)) %>%
    ungroup()

  get_pal <- function(n, pal = "Dark2") {
    nmax <- brewer.pal.info$maxcolors[which(rownames(brewer.pal.info) == pal)]
    if(n <= nmax)
      cols <- brewer.pal(nmax, pal)[seq_len(n)]
    else
      cols <- colorRampPalette(brewer.pal(nmax, pal))(n)
    return(cols)
  }

  ggplot() +
    geom_area(
      data = df1,
      aes(week_one, y = mid, fill = fct_rev(strain)),
      color = 'black'
    ) +
    {if(errorbars)
       geom_errorbar(
         data = df1_adj,
         aes(week_one, ymin = lower, ymax = upper),
         width = 0
       )} +
    {if(errorbars) geom_point(
      data = df1_adj,
      aes(week_one, y = mid),
      size = 1
    )} +
    scale_y_continuous(
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y"
    ) +
    scale_fill_manual(
      values = get_pal(length(variants), palette),
      labels = setNames(get_var_lab(variants, swap = TRUE),
                        paste0("str", seq_len(length(variants))))
    ) +
    scale_color_manual(
      values = rev(get_pal(length(variants), palette)),
      labels = setNames(get_var_lab(variants, swap = TRUE),
                        paste0("str", seq_len(length(variants)))),
      ) +
    labs(x = NULL, y = "Proportion of sequences", fill = "Variant") +
    {if(length(countries) > 1) facet_wrap(~ report_country, ncol = 5)} +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom'
    )

}

## visualise fits
vis_lineage <- function(gis,
                        weeks = 2,
                        log = FALSE,
                        date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                        min_date = as.Date("2021-03-01"),
                        top_n = 30,
                        min_sequences = 500,
                        palette = "Dark2",
                        n_strains = 20) {

  fits <- nest(gis, data = -report_country) %>%
    filter(map_dbl(data, nrow) > min_sequences)

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>%
    arrange(desc(map_dbl(data, nrow))) %>%
    slice(1:top_n)

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    mutate(
      dist_one = warp_distance(date, "week", every = weeks),
      ## pango_lineage = paste0(pango_lineage, " / ", clade),
      pango_lineage = clade
    ) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(as.data.frame(prop.table(table(.$pango_lineage)))) %>%
    rename(strain = Var1, mid = Freq) %>%
    ungroup()

  levs <- df1 %>%
    group_by(report_country, strain) %>%
    filter(week_one >= min_date) %>%
    summarise(mid = sum(mid)) %>%
    arrange(desc(mid)) %>%
    ungroup() %>%
    slice(1:n_strains) %>%
    ungroup() %>%
    mutate(strain = fct_reorder(as.character(strain), mid, mean, .desc = TRUE)) %$%
    c(levels(strain), "Other")

  df1 %<>%
    mutate(
      strain = as.character(strain),
      strain = replace(strain, !strain %in% levs, "Other"),
      strain = factor(strain, levs)
    ) %>%
    group_by(report_country, week_one, strain) %>%
    summarise(mid = sum(mid)) %>%
    complete(report_country, week_one, strain, fill = list(mid = 0))

  get_pal <- function(n, pal = "Dark2") {
    nmax <- brewer.pal.info$maxcolors[which(rownames(brewer.pal.info) == pal)]
    cols <- colorRampPalette(brewer.pal(nmax, pal))(n)
    return(cols)
  }

  pal <- if(length(unique(df1$strain)) <= 8) scale_fill_brewer(palette = palette)
         else scale_fill_manual(values = get_pal(length(unique(df1$strain)), palette))

  ggplot() +
    geom_area(
      data = df1,
      aes(week_one, y = mid, fill = strain),
      color = 'black'
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month"
    ) +
    labs(x = NULL, y = "Proportion of sequences", fill = "Lineage") +
    ## facet_wrap(~ report_country, ncol = 5) +
    pal +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom'
    )

}

## visualise fits
vis_props <- function(fits,
                      countries = unique(fits$report_country),
                      weeks = 2,
                      log = FALSE,
                      date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                      top_n = 30,
                      threshold = 0.01,
                      ncol = 5,
                      days_forecast = 0,
                      font = NULL,
                      base_size = 9,
                      axis_text_size = 6,
                      scaler = 1,
                      rem_most_recent = c("South Africa", "India", "United Kingdom"),
                      split = FALSE) {


  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  ## define bold labels
  txt <- modify_if(
    get_var_lab(fits$labs[[1]], md = TRUE, swap = TRUE),
    ~ grepl(paste0(c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2"), collapse = "|"), get_var_lab(.x)),
    ## ~ .x %in% c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2"),
    ~ paste0("**", .x, "**")
  )

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>%
    filter(!map_lgl(forecast, is.null) & report_country %in% countries)

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      ## arrange(desc(pseudo)) %>%
      slice(1:top_n)
  }

  ## get variant factor levels
  levs <- head(names(fits$labs[[1]]), -1)

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %T>%
    {n_var <<- length(levels(.$y))} %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  ## transform
  df1 %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, levs)) %>%
    filter(mid > 0)

  ## forecast results
  df2 <- unnest(fits, forecast) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    piv_long() %>%
    mutate(strain = factor(strain, levs)) %>%
    filter(mid > 0)

  df1$strain %<>% droplevels()
  df2 %<>%
    filter(as.character(strain) %in% levels(df1$strain)) %>%
    mutate(strain = droplevels(strain))

  ## set palette
  pal <- ifelse(nlevels(df2$strain) > 8, "Set3", "Dark2")

  ## define manual colours
  cols <- get_var_pal(fits$labs[[1]], pal = "Tableau 10", pkg = "tableau")

  if(split) {
    ## df1 %<>% mutate(lower = mid, upper = mid)
    wrap <- facet_grid(
      strain ~ report_country,
      scales = "free",
      labeller = labeller(strain = get_var_lab(fits$labs[[1]])),
      ncol = 5
    )
    guide <- guides(fill = FALSE, color = FALSE)
  } else {
    wrap <- facet_wrap( ~ report_country, ncol = ncol)
    guide <- guides(
      fill = guide_legend(reverse = TRUE),
      color = guide_legend(reverse = TRUE)
    )
  }

  min_date <- df1 %>%
    group_by(report_country) %>%
    summarise(min_date = min(week_one))

  max_date_str <- df1 %>%
    group_by(strain) %>%
    summarise(max_date_str = max(week_one) + days_forecast)

  df2 %<>%
    left_join(min_date, "report_country") %>%
    left_join(max_date_str, "strain") %>%
    filter(date >= min_date & date <= max_date_str)

  df1 %<>%
    group_by(report_country) %>%
    filter(!(week_one == max(week_one) & report_country %in% rem_most_recent))

  ## if(!is.null(font)) {
  ##   font_import()
  ##   loadfonts(device = "win")
  ## }

  ggplot() +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper, group = strain, fill = fct_rev(strain)),
      alpha = 0.5,
      color = NA
    ) +
    geom_line(
      data = df2,
      aes(date, mid, group = strain, color = fct_rev(strain)),
      size = 0.3 * scaler
    ) +
    geom_line(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.3 * scaler,
      alpha = 0.5
    ) +
    geom_point(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.7 * scaler,
      alpha = 1
    ) +
    geom_errorbar(
      data = df1,
      aes(week_one, ymin = lower, ymax = upper, color = fct_rev(strain)),
      width = 0,
      alpha = 0.5,
      size = 0.5 * scaler
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent,
      limits = c(0, 1)
    ) +
    scale_x_date(
      expand = c(0, 0),
      limits = range(df2$date),
      date_breaks = "2 month",
      date_labels = "%b\n%Y"
    ) +
    scale_fill_manual(
      values = cols,
      name = "Variant",
      labels = txt
    ) +
    scale_color_manual(
      values = cols,
      name = "Variant",
      labels = txt
    ) +
    labs(x = NULL, y = "Proportion of sequences") +
    guide +
    wrap +
    theme_minimal(base_size = base_size) +
    theme(
      legend.position = 'bottom',
      legend.text = element_markdown(family = font),
      panel.grid = element_blank(),
      axis.text.x = element_text(size = axis_text_size, family = font),
      text = element_text(family = font)
    )

}

## visualise fits
vis_props_single <- function(fits,
                      countries = unique(fits$report_country),
                      weeks = 2,
                      log = FALSE,
                      date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                      top_n = 30,
                      title = NULL,
                      byrow = TRUE,
                      split = !missing(countries) & length(countries) > 1,
                      threshold = 0.01) {


  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  ## define bold labels
  txt <- modify_if(
    fits$labs[[1]],
    ~ .x %in% c("P.1", "B.1.1.7", "B.1.351"),
    ~ paste0("**", .x, "**")
  )

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>%
    filter(!map_lgl(forecast, is.null) & report_country %in% countries)

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      ## arrange(desc(pseudo)) %>%
      slice(1:top_n)
  }

  ## get variant factor levels
  levs <- head(names(fits$labs[[1]]), -1)

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %T>%
    {n_var <<- length(levels(.$y))} %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  ## transform
  df1 %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, levs)) %>%
    filter(mid > 0)

  ## forecast results
  df2 <- unnest(fits, forecast) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    piv_long() %>%
    mutate(strain = factor(strain, levs)) %>%
    filter(mid > 0)

  df1$strain %<>% droplevels()
  df2 %<>%
    filter(as.character(strain) %in% levels(df1$strain)) %>%
    mutate(strain = droplevels(strain))

  ## set palette
  pal <- ifelse(nlevels(df2$strain) > 8, "Set3", "Dark2")

  ## define manual colours
  cols <- get_var_pal(fits$labs[[1]])

  min_date <- df1 %>%
    group_by(report_country) %>%
    summarise(min_date = min(week_one))

  max_date_str <- df1 %>%
    group_by(strain) %>%
    summarise(max_date_str = max(week_one))

  df2 %<>%
    left_join(min_date, "report_country") %>%
    left_join(max_date_str, "strain") %>%
    filter(date >= min_date & date <= max_date_str)

  if(split) {
    if(byrow) {
      wrap <- facet_grid(report_country ~ strain, labeller = labeller(strain = fits$labs[[1]]))
    } else {
      wrap <- facet_grid(strain ~ report_country, labeller = labeller(strain = fits$labs[[1]]))
    }
    guide <- guides(fill = FALSE, color = FALSE)
  } else {
    wrap <- facet_wrap( ~ strain, nrow = 1, labeller = labeller(strain = fits$labs[[1]]))
    guide <- guides(fill = FALSE, color = FALSE)
  }

  ggplot() +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper, group = strain, fill = fct_rev(strain)),
      alpha = 0.5,
      color = NA
    ) +
    geom_line(
      data = df2,
      aes(date, mid, group = strain, color = fct_rev(strain)),
      size = 0.3
    ) +
    geom_line(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.3,
      alpha = 0.5
    ) +
    geom_point(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.7,
      alpha = 1
    ) +
    geom_errorbar(
      data = df1,
      aes(week_one, ymin = lower, ymax = upper, color = fct_rev(strain)),
      width = 0,
      alpha = 0.5,
      size = 0.5
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0.01, 0),
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.02, 0),
      ## limits = range(df2$date),
      breaks = seq.Date(as.Date("2020-11-01"), as.Date("2021-03-01"), by = "2 months"),
      ## date_breaks = "2 month",
      date_labels = "%b\n%Y"
    ) +
    scale_fill_manual(
      values = cols,
      name = "Variant",
      labels = txt
    ) +
    scale_color_manual(
      values = cols,
      name = "Variant",
      labels = txt
    ) +
    labs(x = NULL, y = "Proportion of sequences", title = title) +
    guide +
    wrap +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = 'bottom',
      legend.text = element_markdown(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6)
    )

}

## visualise fits
vis_country_proportions <- function(data,
                                    labs,
                                    forecast,
                                    weeks = 2,
                                    split = TRUE,
                                    db_date = NULL,
                                    date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                                    title = NULL,
                                    days_forecast = 14
                                    ) {

  ## pivot longer
  piv_long <- function(x) {
    pivot_longer(
      x, contains(c("mid", "lower", "upper")),
      names_to = c(".value", "strain"),
      names_pattern = "(.+)_(.+)"
    )
  }

  ## get variant factor levels
  levs <- head(names(labs), -1)

  ## extract empirical estimates
  df1 <- data %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(get_multinom_ci(.$y)) %>%
    pivot_wider(names_from = "strain", values_from = c("lower", "mid", "upper")) %>%
    ungroup()

  ## transform
  df1 %<>% piv_long() %>% drop_na() %>%
    mutate(strain = factor(strain, levs)) %>%
    filter(mid > 0)

  ## forecast results
  df2 <- forecast %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    piv_long() %>%
    mutate(strain = factor(strain, levs)) %>%
    filter(mid > 0)

  df1$strain %<>% droplevels()
  df2 %<>%
    filter(as.character(strain) %in% levels(df1$strain)) %>%
    mutate(strain = droplevels(strain))

  ## define manual colours
  cols <- get_var_pal(labs, pal = "Set1")

  min_date <- df1 %>%
    group_by(report_country) %>%
    summarise(min_date = min(week_one))

  max_date_str <- df1 %>%
    group_by(strain) %>%
    summarise(max_date_str = max(week_one) + days_forecast)

  df2 %<>%
    mutate(report_country = data$report_country[1]) %>%
    left_join(min_date, "report_country") %>%
    left_join(max_date_str, "strain") %>%
    filter(date >= min_date & date <= max_date_str)

  if(split) {
    facet <- facet_wrap( ~ strain, nrow = 1, labeller = labeller(strain = get_var_lab(labs)))
    guide <- guides(fill = FALSE, color = FALSE)
  } else {
    facet <- guide <- NULL
  }

  ## if(length(unique(df1$week_one)) < 2) return(NULL)

  ggplot() +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper, fill = fct_rev(strain)),
      alpha = 0.5,
      color = NA
    ) +
    geom_line(
      data = df2,
      aes(date, mid, color = fct_rev(strain)),
      size = 0.3
    ) +
    geom_line(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.3,
      alpha = 0.5
    ) +
    geom_point(
      data = df1,
      aes(week_one, y = mid, color = fct_rev(strain)),
      size = 0.7,
      alpha = 1
    ) +
    geom_errorbar(
      data = df1,
      aes(week_one, ymin = lower, ymax = upper, color = fct_rev(strain)),
      width = 0,
      alpha = 0.5,
      size = 0.5
    ) +
    scale_y_continuous(
      expand = c(0.01, 0),
      breaks = seq(0, 1, 0.25),
      limits = c(0, 1),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.02, 0),
      ## breaks = seq.Date(as.Date("2020-10-01"), as.Date("2021-05-01"), by = "2 months"),
      ## date_breaks = "2 month",
      limits = as.Date(c(NA, Sys.Date())),
      date_labels = "%b\n%Y"
    ) +
    scale_fill_manual(
      values = cols,
      name = "Variant",
      labels = get_var_lab(labs, md = TRUE)
    ) +
    scale_color_manual(
      values = cols,
      name = "Variant",
      labels = get_var_lab(labs, md = TRUE)
    ) +
    labs(
      x = NULL,
      y = "Proportion of sequences",
      title = title,
      caption = glue(
        "*Produced by WHO HQ COVID-19 analytics team | Data downloaded from ",
        "GISAID on {format(db_date, '%B %d %Y')}*"
      )
    ) +
    facet +
    guide +
    theme_minimal(base_size = 9) +
    theme(
      ## legend.position = 'bottom',
      legend.text = element_markdown(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6),
      plot.caption = element_markdown(size = 7)
    )

}

## visualise fits for california
vis_glm_place <- function(fits, max_range = 5, min_rsq = 0.02, weeks = 2, r_in_lab = FALSE) {

  fits %<>%
    group_by(strains, place) %>%
    mutate(
      range = map_dbl(tidied, ~ .x$growth_diff[3] - .x$growth_diff[1]),
      est = map_chr(
        tidied,
        function(x) {
          perc <- scales::percent(x$r_ratio_gamma, 1)
          glue("{ifelse(x$r_ratio_gamma[2] > 0, '+', '')}{perc[2]} [{perc[1]} - {perc[3]}]")
        }
      )
    ) %>%
    filter(!is.na(range) & range < max_range & rsq > min_rsq)

  if(r_in_lab) {
    fits %<>% mutate(place = paste(place, est, sep = "\n"))
    text_1 <- text_2 <- NULL
  }

  df1 <- unnest(fits, data)

  df2 <- unnest(fits, forecast)

  df3 <- unnest(fits, data) %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(strains, place, week_one) %>%
    do(binom::binom.exact(sum(.$y), nrow(.))) %>%
    rename(prop_n501 = mean, lower_n501 = lower, upper_n501 = upper)

  if(!r_in_lab) {
    text_1 <- geom_label(
      data = fits,
      aes(label = est),
      x = min(df2$date) + 40,
      y = 0.9,
      color = NA,
      size = 2.5
    )
    text_2 <- geom_text(
      data = fits,
      aes(label = est),
      x = min(df2$date) + 40,
      y = 0.9,
      size = 3
    )
  }

  ggplot() +
    ## geom_vline(xintercept = Sys.Date(), linetype = 2) +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    geom_line(
      data = df2,
      aes(date, y = mid)
    ) +
    geom_point(
      data = df3,
      aes(week_one, y = prop_n501),
      size = 0.8,
      alpha = 1
    ) +
    geom_errorbar(
      data = df3,
      aes(week_one, ymin = lower_n501, ymax = upper_n501),
      width = 0,
      size = 0.4
    ) +
    text_1 +
    text_2 +
    labs(x = NULL, y = NULL) +
    facet_grid(
      strains ~ place,
      switch = "y",
      labeller = labeller(strains = function(x) glue("Proportion {x}"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.placement = "outside"
    )

}

## visualise fits for california
vis_glm <- function(fits, max_range = 5, min_rsq = 0.02, weeks = 2, r_in_lab = TRUE) {

  fits %<>%
    group_by(report_country) %>%
    mutate(
      range = map_dbl(tidied, ~ .x$growth_diff[3] - .x$growth_diff[1]),
      est = map_chr(
        tidied,
        function(x) {
          perc <- scales::percent(x$r_ratio_gamma, 1)
          glue("{ifelse(x$r_ratio_gamma[2] > 0, '+', '')}{perc[2]} [{perc[1]} - {perc[3]}]")
        }
      )
    ) %>%
    filter(!is.na(range) & range < max_range & rsq > min_rsq)

  if(r_in_lab) {
    fits %<>% mutate(place = paste(report_country, est, sep = "\n"))
    text_1 <- text_2 <- NULL
  }

  df1 <- unnest(fits, data)

  df2 <- unnest(fits, forecast)

  df3 <- unnest(fits, data) %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(place, report_country, week_one) %>%
    do(binom::binom.exact(sum(.$y), nrow(.))) %>%
    rename(prop_n501 = mean, lower_n501 = lower, upper_n501 = upper)

  if(!r_in_lab) {
    text_1 <- geom_label(
      data = fits,
      aes(label = est),
      x = min(df2$date) + 40,
      y = 0.9,
      color = NA,
      size = 2.5
    )
    text_2 <- geom_text(
      data = fits,
      aes(label = est),
      x = min(df2$date) + 40,
      y = 0.9,
      size = 3
    )
  }

  levs <- unique(df2$place)
  df1 %<>% mutate(place = factor(place, levs))
  df2 %<>% mutate(place = factor(place, levs))
  df3 %<>% mutate(place = factor(place, levs))

  ggplot() +
    ## geom_vline(xintercept = Sys.Date(), linetype = 2) +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    geom_line(
      data = df2,
      aes(date, y = mid)
    ) +
    geom_point(
      data = df3,
      aes(week_one, y = prop_n501),
      size = 0.8,
      alpha = 1
    ) +
    geom_errorbar(
      data = df3,
      aes(week_one, ymin = lower_n501, ymax = upper_n501),
      width = 0,
      size = 0.4
    ) +
    scale_x_date(date_labels = "%B\n%Y") +
    scale_y_continuous(labels = function(x) scales::percent(x, 1)) +
    text_1 +
    text_2 +
    labs(x = NULL, y = "Proportion VOC") +
    facet_wrap(
      ~ place
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank()
    )

}

## visualise fits
vis_infectiousness <- function(fits, max_range = 5, min_rsq = 0.02, weeks = 2) {

  fits %<>%
    group_by(report_country) %>%
    mutate(
      range = map_dbl(tidied, ~ .x$growth_diff[3] - .x$growth_diff[1]),
      est = map_chr(
        tidied,
        function(x) {
          perc <- scales::percent(x$r_ratio_gamma, 1)
          glue("{perc[2]} [{perc[1]} - {perc[3]}]")
        }
      ),
      lab = glue("{report_country}\n{est}")
    ) %>%
    filter(!is.na(range) & range < max_range & rsq > min_rsq)

  df1 <- unnest(fits, data)

  df2 <- fits %>%
    mutate(
      forecast = map2(
        forecast,
        tidied,
        ~ mutate(
          .x,
          lower = lower*.y$r_ratio_gamma[1],
          mid = mid*.y$r_ratio_gamma[2],
          upper = upper*.y$r_ratio_gamma[3]
        )
      )
    ) %>%
    unnest(forecast) %>%
    ungroup()

  df2 <- tibble(
    report_country = df2$report_country,
    lab = df2$lab,
    date = df2$date,
    bind_rows(as.data.frame(t(apply(select(df2, lower:upper), 1, sort))))
  ) %>%
    setNames(c("report_country", "lab", "date", "lower", "mid", "upper"))

  ggplot() +
    geom_vline(xintercept = Sys.Date(), linetype = 2) +
    geom_ribbon(
      data = df2,
      aes(date, ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    geom_line(
      data = df2,
      aes(date, y = mid)
    ) +
    scale_y_continuous(
      name = expression(Change~R[t]),
      labels = scales::percent
    ) +
    labs(x = NULL) +
    facet_wrap(~ lab) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank()
    )

}

## extract model outputs
extract <- function(model, r_wt = 1,
                    forecast = NULL,
                    type = "nls",
                    gen_time = 5.2, gen_time_sd = 1.9, gen_time_voc = 5.2) {


  if(type == "mglm") {

    out <- tidy_glm(model) %>%
      transmute(
        strain = y.level,
        term = term,
        lower = estimate - 1.96*std.error,
        mid = estimate,
        upper = estimate + 1.96*std.error,
        se = std.error
      ) %>%
      pivot_longer(lower:upper, names_to = "conf") %>%
      pivot_wider(names_from = term, values_from = c(se, value)) %>%
      rename(
        start = `value_(Intercept)`,
        growth_diff = value_t,
        start_se = `se_(Intercept)`,
        growth_diff_se = se_t,
      )

    missing <- setdiff(na.omit(str_aft(names(forecast), "mid_")), out$strain)
    missing <- missing[missing != "str1"]

    if(length(missing) > 0) {
      out %<>%
        bind_rows(
          tibble(
            strain = rep(missing, each = 3),
            conf = rep(c("lower", "mid", "upper"), length(missing))
          )
        )
    }

  } else {

    out <- tidy_glm(model) %>%
      transmute(
        term = term,
        lower = estimate - 1.96*std.error,
        mid = estimate,
        upper = estimate + 1.96*std.error
      ) %>%
      pivot_longer(lower:upper, names_to = "conf") %>%
      pivot_wider(names_from = term) %>%
      rename(start = `(Intercept)`, growth_diff = t)

  }

  out %<>%
    mutate(
      ## calculate start proportion from GLM intercept
      across(c(start, start_se), ~1/(1 + exp(-.x))),
      across(
        c(growth_diff, growth_diff_se),
        list(
          delta = ~ exp(.x*gen_time) - 1,
          exponential = ~ .x*gen_time/r_wt,
          gamma = ~ get_gam_ratio(.x, r_wt = r_wt, w_wt = gen_time, w_voc = gen_time_voc, w_sd = gen_time_sd)
        ),
        .names = "r_ratio_{.fn}{str_remove(.col, 'growth_diff')}"
      )
    )

  return(out)

}

## calculate the ratio of reproduction numbers assuming a gamma distributed
## generation time distribution (see Park et al
## https://www.sciencedirect.com/science/article/pii/S1755436518300847)
get_gam_ratio <- function(growth_diff, r_wt = 1, w_wt = 5.2, w_voc = w_wt, w_sd = 1.9, kappa = NULL) {
  if(is.null(kappa)) kappa <- (w_sd/w_wt)^2
  ## (r_wt^kappa + kappa*w_wt*growth_diff)^(1/kappa)/r_wt - 1
  (1 + (r_wt^kappa - 1)*w_voc/w_wt + kappa*w_voc*growth_diff)^(1/kappa)/r_wt - 1
}

## calculate the growth rate from reproduction number, generation time distribution
calc_growth <- function(r_wt = 1, w_wt = 5.2, w_sd_wt = 1.9, r_diff = 0, w_diff = 0, w_sd_diff = 0, mode = "add") {
  w <- w_wt + w_diff
  w_sd <- w_sd_wt + w_sd_diff
  kappa <- (w_sd/w)^2
  r <- if(mode == "add") r_wt + r_diff else r_wt*(1 + r_diff)
  (r^kappa - 1)/(kappa*w)
}

## estimate r_wt from the proportion voc, growth diff and observed r
est_r_wt <- function(model, country, forecast, index = 1, epinow,
                     gen_time = 5.2, gen_time_sd = 1.9,
                     prop_range = c(0.01, 0.99),
                     w_diff_sq = seq(-0.3, 0.3, length = 50),
                     run_r_wt = TRUE,
                     type = "glm") {

  message(scales::percent(index, 1))

  if(run_r_wt) {

    if(type == "glm") {

      tidied <- tidy(model) %>%
        transmute(
          term = term,
          lower = estimate - 1.96*std.error,
          mid = estimate,
          upper = estimate + 1.96*std.error
        ) %>%
        pivot_longer(lower:upper, names_to = "conf") %>%
        pivot_wider(names_from = term) %>%
        setNames(c("conf", "start", "growth_diff"))

      ## combine empirical rt estimates with prop_voc for times when when prop_voc
      ## is between 5 and 95%
      joint <- filter(forecast, mid > prop_range[1] & mid < prop_range[2]) %>%
        inner_join(
          filter(epinow, report_country == country, type != "forecast"),
          "date"
        ) %>%
        transmute(
          date, lower, mid, upper,
          lower_r_tot = lower_90,
          mid_r_tot = median,
          upper_r_tot = upper_90
        )

      ## estimate r_wt from the proportion voc, growth diff and observed r
      infer_r <- function(prop_voc, r_tot, growth_diff, r_wt_sq = seq(0.01, 4, length = 1000)) {
        ## this is solved for out == 0
        out <- prop_voc*r_wt_sq *
          (1 + get_gam_ratio(growth_diff, r_wt = r_wt_sq, w_sd = gen_time_sd)) +
          (1 - prop_voc)*r_wt_sq - r_tot
        ## return estimate
        ind <- which.min(abs(out))
        if(ind %in% c(1, length(out))) {
          stop("Rt estimate has not converged")
        }
        return(r_wt_sq[ind])
      }

      pmap(
        list(
          select(joint, lower:upper),
          tidied$growth_diff,
          select(joint, contains("r_tot"))
        ),
        function(prop_voc, growth_diff, r_tot) {
          map2_dfr(
            prop_voc, r_tot,
            ~ tibble(
              prop_voc = .x,
              r_wt = infer_r(.x, .y, growth_diff = growth_diff),
              r_voc = r_wt*(1 + get_gam_ratio(growth_diff,
                                              r_wt = r_wt,
                                              w_sd = gen_time_sd)),
              ratio = r_voc/r_wt - 1
            )
          )
        }
      ) %>%
        data.frame() %>%
        rename_all(str_replace, "\\.", "_") %>%
        mutate(
          date = joint$date,
          mid_r_tot = joint$mid_r_tot,
          lower_r_tot = joint$lower_r_tot,
          upper_r_tot = joint$upper_r_tot
        ) %>%
        select(date, everything()) %>%
        return()

    } else if(type == "mglm") {

      tidied <- tidy_glm(model) %>%
        transmute(
          strain = y.level,
          term = term,
          lower = estimate - 1.96*std.error,
          mid = estimate,
          upper = estimate + 1.96*std.error
        ) %>%
        pivot_longer(lower:upper, names_to = "conf") %>%
        pivot_wider(names_from = term) %>%
        rename(start = `(Intercept)`, growth_diff = t)

      missing <- setdiff(na.omit(str_aft(names(forecast), "mid_")), tidied$strain)
      missing <- missing[missing != "str1"]

      if(length(missing) > 0) {
        tidied %<>%
          bind_rows(
            tibble(
              strain = rep(missing, each = 3),
              conf = rep(c("lower", "mid", "upper"), length(missing)),
              start = NA,
              growth_diff = NA
            )
          )
      }

      ## combine empirical rt estimates with prop_voc for times when when prop_voc
      ## is between 5 and 95%

      epi <- filter(epinow, report_country == country, type != "forecast")
      joint <- filter(forecast, mid_str1 < (1 - prop_range[1])) %>%
        inner_join(epi, "date") %>%
        select(date, contains("str"), lower_90, median, upper_90) %>%
        rename(
          lower_r_tot = lower_90,
          mid_r_tot = median,
          upper_r_tot = upper_90
        )

      if(nrow(joint) == 0) stop("No Rt estimates available")

      ## estimate r_wt from the proportion voc, growth diff and observed r
      infer_r_old <- function(prop_voc, r_tot, growth_diff,
                          w_diff = seq(-0.3, 0.3, length = 25)
                          ) {

        ## make sure prop_voc sums to 1
        adj <- 1 - sum(prop_voc[1,-1], na.rm = TRUE)
        prop_voc[1,1] <- adj
        if(adj < 0) {
          prop_voc[1,1] <- 0
          prop_voc <- prop_voc/sum(prop_voc)
        }

        ## add a growth_diff of 0 for str1 (i.e. WT)
        growth_diff <- bind_rows(
          tibble(strain = "str1", conf = NA, start = NA, growth_diff = 0),
          growth_diff
        )

        ## ## calculate the mean Rt given WT rt and proportion of VOCs
        ## mean_r <- map2(
        ##   growth_diff$growth_diff,
        ##   prop_voc,
        ##   ~ .y * r_wt_sq * (1 + get_gam_ratio(.x, r_wt = r_wt_sq, w_sd = gen_time_sd))
        ## ) %>%
        ##   pmap_dbl(function(...) sum(..., na.rm = TRUE))

        ## increase Rt parameter space search with time
        error <- TRUE
        range <- 0.6
        step <- 0
        while(error) {

          ## define r_wt search space
          sq <- 1 + seq(range*step, range*(step + 1), length = 20)
          r_wt_sq <- c(rev(1/sq), sq)*r_tot

          ## define combination of r_wt and w_voc
          vars <- expand.grid(r_wt = r_wt_sq, w_voc = gen_time*(1 + w_diff))

          ## calculate the mean Rt given WT rt and proportion of VOCs
          mean_r <- map2(
            growth_diff$growth_diff,
            prop_voc,
            function(growth_diff, prop_voc) {
              pmap_dbl(
                vars,
                function(r_wt, w_voc) {
                  prop_voc * r_wt * (1 + get_gam_ratio(growth_diff, r_wt = r_wt, w_voc = w_voc, w_sd = gen_time_sd))
                }
              )
            }
          ) %>%
            pmap_dbl(function(...) sum(..., na.rm = TRUE))

          ## ## calculate the mean Rt given WT rt and proportion of VOCs
          ## ratio <- map2(
          ##   growth_diff$growth_diff,
          ##   prop_voc,
          ##   function(growth_diff, prop_voc) {
          ##     pmap_dbl(
          ##       vars,
          ##       function(r_wt, w_voc) {
          ##         get_gam_ratio(growth_diff, r_wt = r_wt, w_voc = w_voc, w_sd = gen_time_sd)
          ##       }
          ##     )
          ##   }
          ## )

          ## get the best guess for r_wt (i.e. calculated r_tot is closest to observed value)
          best <- vars %>%
            mutate(est_diff = abs(mean_r - r_tot)) %>%
            left_join(tibble(w_voc = gen_time*(1 + w_diff), w_diff = w_diff), "w_voc") %>%
            group_by(w_diff) %>%
            summarise(r_wt = r_wt[which.min(est_diff)])

          error <- any(best$r_wt %in% r_wt_sq[c(1, length(r_wt_sq))])
          step <- step + 1

          if(step > 100) {
            error <- FALSE
            message("Rt estimates didn't converge")
          }

        }

        ## return best estimates
        return(best)

      }

      calculate_r <- function(r_wt, r_tot, w_diff, growth_diff, prop_voc) {

        map2(
          growth_diff,
          prop_voc,
          function(growth_diff, prop_voc) {
            map_dbl(
              r_wt,
              function(r_wt) {
                prop_voc * r_wt * (1 + get_gam_ratio(growth_diff, r_wt = r_wt, w_voc = gen_time*(1 + w_diff)))
              }
            )
          }
        ) %>%
          pmap_dbl(function(...) sum(..., na.rm = TRUE)) - r_tot

      }

      infer_r <- function(prop_voc, r_tot, growth_diff, w_diff = seq(-0.3, 0.3, length = 25)) {

        ## make sure prop_voc sums to 1
        adj <- 1 - sum(prop_voc[1,-1], na.rm = TRUE)
        prop_voc[1,1] <- adj
        if(adj < 0) {
          prop_voc[1,1] <- 0
          prop_voc <- prop_voc/sum(prop_voc)
        }

        ## add a growth_diff of 0 for str1 (i.e. WT)
        growth_diff <- bind_rows(
          tibble(strain = "str1", conf = NA, start = NA, growth_diff = 0),
          growth_diff
        )

        tibble(
          w_diff = w_diff,
          r_wt = map_dbl(
            w_diff,
            ~ uniroot(
              calculate_r,
              interval = c(0, 2),
              extendInt = "yes",
              r_tot = r_tot,
              w_diff = .x,
              growth_diff = growth_diff$growth_diff,
              prop_voc = unlist(prop_voc)
            )$root
          )
        )

      }

      ## calculate Rt across lower mid and upper growth rate estiamtes, across a
      ## range of generation time differences
      pmap(
        list(
          ## PROP VOC is REVERSED because a lower prop VOC gives a HIGHER R_wt
          ## growth_diff is also reversed because a higher growth difference for
          ## VOC means they have a lower growth diff and r_tot
          map(c("upper_str", "mid_str", "lower_str"), ~ select(joint, contains(.x))),
          map(c("upper", "mid", "lower"), ~ filter(tidied, conf == .x)),
          select(joint, contains("r_tot"))
        ),
        function(prop_voc, growth_diff, r_tot) {
          map2_dfr(
            split(prop_voc, seq_len(nrow(prop_voc))), r_tot,
            function(prop_voc, r_tot) {
              r_str1 <- infer_r(prop_voc, r_tot, growth_diff = growth_diff, w_diff = w_diff_sq)
              df <- data.frame(
                w_diff = r_str1$w_diff,
                prop = prop_voc,
                map_dfc(
                  setNames(c(0, growth_diff$growth_diff), c("str1", growth_diff$strain)),
                  ~ r_str1$r_wt*(1 + get_gam_ratio(.x, r_str1$r_wt,
                                                   w_voc = gen_time*(1 + r_str1$w_diff),
                                                   w_sd = gen_time_sd))
                ) %>%
                rename_all(~ paste0("r_", str_bef(names(prop_voc)[1], "_"), "_", .x))
              ) %>%
                rename_all(str_replace_all, "\\.", "_")
            }
          )
        }
      ) %>%
        data.frame() %>%
        mutate(
          date = rep(joint$date, each = length(w_diff_sq)),
          r_lower_tot = rep(joint$lower_r_tot, each = length(w_diff_sq)),
          r_mid_tot = rep(joint$mid_r_tot, each = length(w_diff_sq)),
          r_upper_tot = rep(joint$upper_r_tot, each = length(w_diff_sq))
        ) %>%
        select(date, everything(), -contains("w_diff.")) %>%
        return()

    }

  } else return(NULL)

}

## hacky rsquared for nls
get_nls_rsq <- function(model, data) {
  f <- predict(model)
  y <- data$y
  1 - sum((y - f)^2)/(length(y)*var(y))
}

## get Rsq equivalent for glm (https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r)
get_glm_rsq <- function(model) {
  with(summary(model), 1 - deviance/null.deviance)
}

## histogram of R ratio estimates
vis_hist <- function(fits, min_rsq = 0.02, max_range = 5, facet = FALSE) {

  fits %<>%
    mutate(range = map_dbl(tidied, ~ .x$growth_diff[3] - .x$growth_diff[1])) %>%
    filter(rsq > min_rsq & range < max_range) %>%
    unnest(tidied) %>%
    ungroup() %>%
    select(report_country, conf, contains("r_ratio")) %>%
    pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
    pivot_wider(names_from = "conf", values_from = "value") %>%
    mutate(
      report_country = fct_reorder(report_country, mid),
      r_model = factor(r_model, c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta"))
    )

  if(!facet) {
    fits %<>% filter(r_model == "r_ratio_gamma")
    color_scale <- scale_color_manual(
      values = "black",
      guide = FALSE
    )
  } else {
    color_scale <- scale_color_brewer(
      palette = "Dark2",
      name = "Generation time distribution",
      labels = function(x) str_to_title(str_remove(x, "r_ratio_"))
    )
  }

  fits %>%
    ggplot(aes(report_country, mid)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, color = r_model),
      position = position_dodge(width = 0.5),
      width = 0.5
    ) +
    geom_point(
      aes(color = r_model),
      position = position_dodge(width = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = 2) +
    color_scale +
    scale_y_continuous(labels = scales::percent) +
    labs(
      x = NULL,
      y = "Change in R"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      legend.position = 'bottom',
      legend.title = element_text(vjust = 0.75)
    )

}

## histogram of R ratio estimates
vis_bar_glm <- function(fits, var = "r_ratio", facet = FALSE, lims = TRUE) {

  labs <- fits$labs[[1]]

  if(var == "r_ratio") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, contains("r_ratio")) %>%
      drop_na() %>%
      pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
      pivot_wider(names_from = "conf", values_from = "value") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, mean),
        r_model = factor(
          r_model,
          c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta")
        ),
        strain = fct_reorder(strain, mid, mean),
        strain_color = factor(strain, head(names(labs), -1))
      ) %>%
      filter(grepl("gamma", r_model))

  } else if(var == "growth_diff") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, growth_diff) %>%
      drop_na() %>%
      pivot_wider(names_from = "conf", values_from = "growth_diff") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, mean),
        strain = fct_reorder(strain, mid, mean),
        strain_color = factor(strain, head(names(labs), -1))
      )

  }

  lims <- if(lims) range(fits$mid) else c(NA, NA)

  fits %>%
    ggplot(aes(report_country, mid)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, color = strain_color),
      position = position_dodge(width = 0.5),
      width = 0
    ) +
    geom_point(
      aes(color = strain_color),
      position = position_dodge(width = 0.5)
    ) +
    scale_color_brewer(
      palette = "Dark2",
      name = "Variant",
      labels = labs
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = lims) +
    labs(
      x = NULL,
      y = glue("Change in {ifelse(var == 'r_ratio', 'R', 'growth rate')}")
    ) +
    ## facet_wrap(~ strain, labeller = labeller(strain = labs)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.5),
      legend.position = 'bottom',
      legend.title = element_text(vjust = 0.75)
    )

}

## histogram of R ratio estimates
vis_violin_glm <- function(fits, var = "r_ratio", facet = FALSE, lims = NULL) {

  labs <- fits$labs[[1]]

  if(var == "r_ratio") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, contains("r_ratio")) %>%
      drop_na() %>%
      pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
      pivot_wider(names_from = "conf", values_from = "value") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, mean),
        r_model = factor(
          r_model,
          c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta")
        ),
        strain = fct_reorder(strain, mid, median),
        strain_color = factor(strain, head(names(labs), -1))
      ) %>%
      filter(grepl("gamma", r_model))

  } else if(var == "growth_diff") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, growth_diff) %>%
      drop_na() %>%
      pivot_wider(names_from = "conf", values_from = "growth_diff") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, median),
        strain = fct_reorder(strain, mid, median),
        strain_color = factor(strain, head(names(labs), -1))
      )

  }

  if(is.null(lims)) lims <- range(fits$mid)

  fits %>%
    ggplot(aes(strain, mid)) +
    ## geom_violin(
    ##   aes(strain, mid, fill = strain_color),
    ##   ## draw_quantiles = 0.5,
    ##   bw = 0.05
    ## ) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, group = report_country),
      position = position_dodge(width = 0.5),
      alpha = 0.5,
      color = "darkgrey",
      width = 0
    ) +
    geom_point(
      aes(group = report_country),
      position = position_dodge(width = 0.5),
      size = 0.5,
      color = "darkgrey"
    ) +
    geom_boxplot(
      aes(fill = strain_color),
      width = 0.2,
      size = 0.6,
      alpha = 0.75,
      color = "black",
      outlier.shape = NA
    ) +
    scale_fill_brewer(
      palette = "Set3",
      guide = FALSE,
      labels = labs,
      drop = FALSE
    ) +
    scale_color_brewer(
      palette = "Set3",
      guide = FALSE,
      labels = labs,
      drop = FALSE
    ) +
    scale_x_discrete(labels = labs, drop = FALSE) +
    scale_y_continuous(labels = scales::percent) +
    coord_cartesian(ylim = lims) +
    labs(
      x = "Variant",
      y = glue("Change in {ifelse(var == 'r_ratio', 'R', 'growth rate')}")
    ) +
    ## facet_wrap(~ strain, labeller = labeller(strain = labs)) +
    theme_minimal() +
    theme(
      legend.position = 'bottom',
      legend.title = element_text(vjust = 0.75)
    )

}

## histogram of R ratio estimates
vis_boxplot_glm <- function(fits, meta, var = "r_ratio", facet = FALSE, lims = NULL,
                            base_size = 9, axis_text_size = NULL, scaler = 1) {

  labs <- get_var_lab(fits$labs[[1]], swap = FALSE, remove = TRUE) %>%
    str_replace_all("B.1.427/B.1.429", "B.1.427/\nB.1.429") %>%
    setNames(names(fits$labs[[1]]))

  cols <- get_var_pal(fits$labs[[1]], pal = "Tableau 10", pkg = "tableau")

  if(var == "r_ratio") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, contains("r_ratio")) %>%
      drop_na() %>%
      pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
      pivot_wider(names_from = "conf", values_from = "value") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, mean),
        r_model = factor(
          r_model,
          c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta")
        ),
        strain = fct_reorder(strain, mid, median),
        strain_color = factor(strain, head(names(labs), -1))
      ) %>%
      filter(grepl("gamma", r_model))

    var2 <- paste0(var, "_gamma_sm")
    meta %<>%
      unnest({{var2}}) %>%
      mutate(
        strain = fct_reorder(strain, mid),
        strain_color = factor(strain, head(names(labs), -1))
      )

  } else if(var == "growth_diff") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, growth_diff) %>%
      drop_na() %>%
      pivot_wider(names_from = "conf", values_from = "growth_diff") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, median),
        strain = fct_reorder(strain, mid, median),
        strain_color = factor(strain, head(names(labs), -1))
      )

    var2 <- paste0(var, "_sm")
    meta %<>%
      unnest({{var2}}) %>%
      mutate(
        strain = fct_reorder(strain, mid),
        strain_color = factor(strain, head(names(labs), -1))
      )

  }

  if(is.null(lims)) lims <- range(fits$mid)

  meta %<>%
    mutate(
      ## voc = fct_rev(ifelse(labs[as.character(strain)] %in%
      ##                      c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2"),
      ##                      "VOC", "VOI")),
      voc = fct_rev(
        ifelse(
          grepl(
            paste0(
              c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2", "Alpha", "Beta", "Gamma", "Delta"),
              collapse = "|"
            ),
            labs[as.character(strain)]
          ),
          "VOC", "VOI"
        )
      )
      ## across(c(mid, lower, upper), function(x) modify_at(x, 1, ~ .x + 0.03))
    )

  fits %>%
    filter(!(strain == "str2" & mid < 0)) %>%
    mutate(
      strain = factor(strain, levels(meta$strain)),
      voc = fct_rev(
        ifelse(
          grepl(
            paste0(
              c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2", "Alpha", "Beta", "Gamma", "Delta"),
              collapse = "|"
            ),
            labs[as.character(strain)]
          ),
          "VOC", "VOI"
        )
      )
    ) %>%
    ggplot(aes(strain, mid)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, group = report_country),
      position = position_dodge(width = 0.5),
      alpha = 1,
      color = "darkgrey",
      size = 0.25 * scaler,
      width = 0
    ) +
    geom_point(
      aes(group = report_country),
      position = position_dodge(width = 0.5),
      size = 0.2 * scaler,
      color = "darkgrey"
    ) +
    geom_errorbar(
      data = meta,
      aes(y = mid, ymin = lower, ymax = upper),
      size = 1 * scaler,
      width = 0
    ) +
    geom_point(
      data = meta,
      aes(y = mid, fill = strain),
      size = 3 * scaler,
      shape = 23
    ) +
    scale_fill_manual(
      values = cols,
      guide = FALSE
    ) +
    scale_color_manual(
      values = cols,
      guide = FALSE
    ) +
    scale_x_discrete(labels = labs) +
    scale_y_continuous(labels = pminus_percent, limits = lims) +
    labs(
      x = NULL,
      y = glue("Change in {ifelse(var == 'r_ratio', 'R', 'growth rate')}")
      ## caption = "* indicates potential immune evasion, ^ indicates increased mortality"
    ) +
    facet_wrap(~ voc, nrow = 1, scales = "free_x") +
    theme_minimal(base_size = base_size) +
    theme(
      strip.background = element_rect(),
      axis.text.x = element_text(size = axis_text_size)
    )

}

## histogram of R ratio estimates
vis_boxplot_sensitivity <- function(fits, meta, var = "r_ratio", facet = FALSE, lims = NULL,
                                    base_size = 9, axis_text_size = NULL, scaler = 1) {

  labs <- get_var_lab(fits$labs[[1]], swap = FALSE, remove = TRUE) %>%
    str_replace_all("B.1.427/B.1.429", "B.1.427/\nB.1.429") %>%
    setNames(names(fits$labs[[1]]))

  cols <- get_var_pal(fits$labs[[1]], pal = "Tableau 10", pkg = "tableau")

  ## fits2 <- fits %>%
  ##   mutate(
  ##     tidied = pmap(
  ##       list(
  ##         glm,
  ##         map(r_wt, possibly(
  ##           ~ filter(.x, abs(w_diff) == min(abs(w_diff))) %$%
  ##             mean(r_mid_str1, na.rm = TRUE),
  ##           1
  ##         )),
  ##         forecast
  ##       ),
  ##       possibly(extract, NULL),
  ##       type = "mglm"
  ##     )
  ##   )

  ## r1 <- map_dbl(fits$r_wt, possibly(
  ##   ~ filter(.x, abs(w_diff) == min(abs(w_diff))) %$%
  ##     mean(r_mid_str1, na.rm = TRUE),
  ##   1
  ## ))

  ## r2 <- map_dbl(fits$r_wt, possibly(
  ##   ~ filter(.x, abs(w_diff + 0.2387755) == min(abs(w_diff + 0.2387755))) %$%
  ##     mean(r_mid_str1, na.rm = TRUE),
  ##   1
  ## ))

  fits1 <- fits %>%
    unnest(tidied) %>%
    ungroup() %>%
    select(report_country, strain, conf, contains("r_ratio")) %>%
    drop_na() %>%
    pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
    pivot_wider(names_from = "conf", values_from = "value") %>%
    mutate(
      report_country = fct_reorder(report_country, mid, mean),
      r_model = factor(
        r_model,
        c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta")
      ),
      strain = fct_reorder(strain, mid, median),
      strain_color = factor(strain, head(names(labs), -1))
    ) %>%
    filter(grepl("gamma", r_model))

  fits2 <- fits1 %>%
    mutate(across(c(lower, mid, upper), ~ .x*0.688))

  ## meta2 <- get_meta(filter(fits2, !report_country %in% c("Turkey", "Zambia")))

  ## fits2 <- fits2 %<>%
  ##   unnest(tidied) %>%
  ##   ungroup() %>%
  ##   select(report_country, strain, conf, contains("r_ratio")) %>%
  ##   drop_na() %>%
  ##   pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
  ##   pivot_wider(names_from = "conf", values_from = "value") %>%
  ##   mutate(
  ##     report_country = fct_reorder(report_country, mid, mean),
  ##     r_model = factor(
  ##       r_model,
  ##       c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta")
  ##     ),
  ##     strain = fct_reorder(strain, mid, median),
  ##     strain_color = factor(strain, head(names(labs), -1))
  ##   ) %>%
  ##   filter(grepl("gamma", r_model)) %>%
  ##   mutate(across(c(lower, mid, upper), ~ x*0.688))

  var2 <- paste0(var, "_gamma_sm")

  meta %<>%
    unnest({{var2}}) %>%
    mutate(
      strain = fct_reorder(strain, mid),
      strain_color = factor(strain, head(names(labs), -1))
    )

  if(is.null(lims)) lims <- range(fits$mid)

  meta %<>%
    mutate(
      voc = fct_rev(
        ifelse(
          grepl(
            paste0(
              c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2", "Alpha", "Beta", "Gamma", "Delta"),
              collapse = "|"
            ),
            labs[as.character(strain)]
          ),
          "VOC", "VOI"
        )
      )
    )

  meta2 <- meta %>% mutate(across(c(lower, mid, upper), ~ .x*0.688))

  ## meta2 %<>%
  ##   mutate(
  ##     voc = fct_rev(
  ##       ifelse(
  ##         grepl(
  ##           paste0(
  ##             c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2", "Alpha", "Beta", "Gamma", "Delta"),
  ##             collapse = "|"
  ##           ),
  ##           labs[as.character(strain)]
  ##         ),
  ##         "VOC", "VOI"
  ##       )
  ##     ),
  ##     across(c(lower, mid, upper), ~ x*0.688)
  ##   )

  meta_join <- bind_rows(
    mutate(meta, mod = "Unchanged serial interval (5.2 days)"),
    mutate(meta2, mod = "Shorter serial interval (4.0 days)")
  ) %>%
    mutate(mod = fct_rev(mod)) %>%
    filter(strain == names(fits$labs[[1]])[which(fits$labs[[1]] == "B.1.617.2")])

  bind_rows(
    mutate(fits1, mod = "Unchanged serial interval (5.2 days)"),
    mutate(fits2, mod = "Shorter serial interval (4.0 days)")
  ) %>%
    filter(strain == names(fits$labs[[1]])[which(fits$labs[[1]] == "B.1.617.2")]) %>%
    mutate(
      strain = factor(strain, levels(meta$strain)),
      mod = fct_rev(mod)
    ) %>%
    ggplot(aes(mod, mid)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, group = report_country),
      position = position_dodge(width = 0.5),
      alpha = 1,
      color = "darkgrey",
      size = 0.35 * scaler,
      width = 0
    ) +
    geom_point(
      aes(group = report_country),
      position = position_dodge(width = 0.5),
      size = 0.4 * scaler,
      color = "darkgrey"
    ) +
    geom_errorbar(
      data = meta_join,
      aes(y = mid, ymin = lower, ymax = upper),
      size = 1.4 * scaler,
      width = 0
    ) +
    geom_point(
      data = meta_join,
      aes(y = mid, fill = strain),
      size = 3 * scaler,
      shape = 23
    ) +
    scale_fill_manual(
      values = cols,
      guide = FALSE
    ) +
    scale_color_manual(
      values = cols,
      guide = FALSE
    ) +
    scale_x_discrete(labels = labs) +
    scale_y_continuous(labels = pminus_percent, limits = lims) +
    labs(
      x = NULL,
      y = glue("Change in {ifelse(var == 'r_ratio', 'R', 'growth rate')}")
      ## caption = "* indicates potential immune evasion, ^ indicates increased mortality"
    ) +
    ## facet_wrap(~ mod, nrow = 1, scales = "free_x") +
    theme_minimal(base_size = base_size) +
    theme(
      strip.background = element_rect(),
      axis.text.x = element_text(size = axis_text_size)
    )

}

## histogram of R ratio estimates
vis_boxplot_single <- function(fits, meta, var = "r_ratio", facet = FALSE, lims = NULL) {

  labs <- fits$labs[[1]]

  cols <- get_var_pal(fits$labs[[1]])

  if(var == "r_ratio") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, contains("r_ratio")) %>%
      drop_na() %>%
      pivot_longer(contains("r_ratio"), names_to = "r_model") %>%
      pivot_wider(names_from = "conf", values_from = "value") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, mean),
        r_model = factor(
          r_model,
          c("r_ratio_exponential", "r_ratio_gamma", "r_ratio_delta")
        ),
        strain = fct_reorder(strain, mid, median),
        strain_color = factor(strain, head(names(labs), -1))
      ) %>%
      filter(grepl("gamma", r_model))

    var2 <- paste0(var, "_gamma_sm")
    meta %<>%
      unnest({{var2}}) %>%
      mutate(
        strain = fct_reorder(strain, mid),
        strain_color = factor(strain, head(names(labs), -1))
      )

  } else if(var == "growth_diff") {

    fits %<>%
      unnest(tidied) %>%
      ungroup() %>%
      select(report_country, strain, conf, growth_diff) %>%
      drop_na() %>%
      pivot_wider(names_from = "conf", values_from = "growth_diff") %>%
      mutate(
        report_country = fct_reorder(report_country, mid, median),
        strain = fct_reorder(strain, mid, median),
        strain_color = factor(strain, head(names(labs), -1))
      )

    var2 <- paste0(var, "_sm")
    meta %<>%
      unnest({{var2}}) %>%
      mutate(
        strain = fct_reorder(strain, mid),
        strain_color = factor(strain, head(names(labs), -1))
      )

  }

  if(is.null(lims)) lims <- range(fits$mid)

  meta %<>% mutate(
    voc = fct_rev(ifelse(labs[as.character(strain)] %in% c("P.1", "B.1.1.7", "B.1.351"),
                         "VOC", "VOI"))
  )

  fits %>%
    mutate(
      strain = factor(strain, levels(meta$strain)),
      voc = fct_rev(ifelse(labs[as.character(strain)] %in% c("P.1", "B.1.1.7", "B.1.351"),
                           "VOC", "VOI"))
    ) %>%
    ggplot(aes(strain, mid)) +
    ## geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, color = strain),
      size = 1,
      width = 0
    ) +
    geom_point(
      aes(group = report_country, color = strain),
      size = 3,
      shape = 18
    ) +
    scale_fill_manual(
      values = cols,
      guide = FALSE
    ) +
    scale_color_manual(
      values = cols,
      guide = FALSE
    ) +
    scale_x_discrete(labels = labs) +
    scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
    labs(
      x = NULL,
      y = glue("Change in {ifelse(var == 'r_ratio', 'R', 'growth rate')}")
    ) +
    facet_wrap(~ voc, nrow = 1, scales = "free_x") +
    theme_minimal(base_size = 9) +
    theme(
      strip.background = element_rect(),
      axis.text.x = element_text(size = 6)
    )

}

## fit models, extract summaries, forecast
get_nls_fits <- function(df) {

  df %>%
    drop_na(week, prop_n501) %>%
    group_by(report_country) %>%
    arrange(week) %>%
    slice(min(which(prop_n501 > 0)):nrow(.)) %>%
    mutate(t = as.numeric(week - min(week)), y = prop_n501) %>%
    nest() %>%
    mutate(
      ## fits models
      growth = map(data, possibly(fit_growth, NULL)),
      r = map(data, possibly(fit_r, NULL)),
      across(
        c(growth, r),
        map,
        possibly(predict, NULL),
        .names = "{.col}_predict"
      ),
      ## summarise estimates
      across(
        c(growth, r),
        map,
        possibly(extract, tibble(term = NA, lower = NA, mid = NA, upper = NA)),
        .names = "{.col}_tidied"
      ),
      ## calculate rsq
      across(
        c(growth, r),
        ~ map2_dbl(.x, data, possibly(get_nls_rsq, NA)),
        .names = "{.col}_rsq"
      ),
      growth_forecast = map2(
        growth_tidied,
        data,
        possibly(forecast_growth, tibble(date = NA, lower = NA, mid = NA, upper = NA))
       ),
      r_forecast = map2(
        r_tidied,
        data,
        possibly(forecast_r, tibble(date = NA, lower = NA, mid = NA, upper = NA)),
        gen_time = 5.0
      )
    )

}

## fit models, extract summaries, forecast
get_glm_fits <- function(gis, epinow, variants = "B.1.1.7",
                         gen_time = 5.2, gen_time_sd = 1.9,
                         t_start = as.Date("2020-09-01"), t_end = Sys.Date() + 45,
                         prop_range = c(0, 1)) {

  gis %>%
    mutate(y = pango_lineage %in% variants) %>%
    group_by(report_country) %>%
    filter(
      sum(y) > 50,
      sum(date > as.Date("2020-11-01")) > 250
    ) %>%
    arrange(date) %>%
    slice(min(which(y)):nrow(.)) %>%
    mutate(t = as.numeric(date - min(date))) %>%
    nest() %>%
    mutate(
      ## fit models
      glm = map(data, possibly(fit_glm, NULL)),
      ## predict from these models
      predict = map(glm, possibly(predict, NULL)),
      ## calculate rsq
      rsq = map_dbl(glm, get_glm_rsq),
      ## forecast proportions
      forecast = map2(glm, data, forecast_glm, t_start = t_start, t_end = t_end),
      ## estimate r_wt
      r_wt = pmap(
        list(glm, report_country, forecast),
        possibly(est_r_wt, NULL),
        epinow = epinow,
        gen_time = gen_time,
        gen_time_sd = gen_time_sd,
        prop_range = prop_range
      ),
      ## extract summary statistics
      tidied = map2(
        glm,
        map(r_wt, ~ mean(.x$mid_r_wt, na.rm = TRUE)),
        possibly(extract, tibble(term = NA, lower = NA, mid = NA, upper = NA)),
        type = "glm"
      )
    ) %>%
    ungroup()

}

## fit models, extract summaries, forecast
get_mglm_fits <- function(gis, epinow, variants = c("B.1.1.7", "B.1.351", "P.1"),
                         gen_time = 5.2, gen_time_sd = 1.9,
                         t_start = as.Date("2020-09-01"), t_end = Sys.Date() + 45,
                         min_n_var = 25, min_n_tot = 250, min_n_since = as.Date("2020-11-01"),
                         run_r_wt = TRUE,
                         save = TRUE,
                         w_diff_sq = seq(-0.3, 0.3, length = 50),
                         prop_range = c(0.01, 0.99)) {

  fits <- gis %>%
    group_by(report_country, pango_lineage) %>%
    mutate(pango_lineage = ifelse(n() > min_n_var, pango_lineage, "Other")) %>%
    ungroup() %>%
    mutate(
      y = factor(
        paste0("str", replace_na(match(pango_lineage, variants) + 1, 1)),
        levels = paste0("str", seq_len(length(variants) + 1))
      )
    ) %T>%
    {nms <<- setNames(c("Other", variants, "Joint"), c(levels(.$y), "tot"))} %>%
    group_by(report_country) %>%
    filter(
      sum(as.numeric(y) > 1) >= min_n_var,
      sum(date > min_n_since) >= min_n_tot
    ) %>%
    arrange(date) %>%
    slice(min(which(as.numeric(y) > 1)):nrow(.)) %>%
    group_by(report_country) %>%
    mutate(t = as.numeric(date - min(date))) %>%
    nest() %>%
    ungroup() %T>%
    {n_countries <<- nrow(.)} %>%
    mutate(
      ## fit models
      glm = map(data, possibly(fit_mglm, NULL)),
      ## predict from these models
      predict = map(glm, possibly(predict, NULL)),
      ## ## calculate rsq
      ## rsq = map_dbl(glm, get_glm_rsq),
      ## forecast proportions
      forecast = map2(glm, data, possibly(forecast_mglm, NULL), t_start = t_start, t_end = t_end),
      ## estimate r_wt
      r_wt = pmap(
        list(glm, report_country, forecast, seq_len(n_countries)/n_countries),
        possibly(est_r_wt, NULL),
        epinow = epinow,
        gen_time = gen_time,
        gen_time_sd = gen_time_sd,
        prop_range = prop_range,
        type = "mglm",
        run_r_wt = run_r_wt,
        w_diff_sq = w_diff_sq
      ),
      ## extract summary statistics
      tidied = pmap(
        list(
          glm,
          if(exists("r_wt"))
            map(r_wt, possibly(
              ~ filter(.x, w_diff == min(abs(w_diff)), r_mid_str1 > (1 - prop_range[2])) %$%
                mean(r_mid_str1, na.rm = TRUE),
              1
            ))
          else rep(1, n_countries),
          forecast
        ),
        possibly(extract, NULL),
        type = "mglm"
      ),
      ## ## ## ## calculate a range of r values given the generation time
      r_range = map2(
        r_wt, tidied,
        possibly(get_r_range, NULL),
        gen_time = gen_time,
        gen_time_sd = gen_time_sd,
        run_r_wt = run_r_wt
      ),
      ## get pseudo R
      pseudo = map_dbl(glm, possibly(DescTools::PseudoR2, NA), "Nagelkerke"),
      ## get strain labels
      labs = rep(list(nms), length(glm)),
      ## get most common variant on the last day of data
      most_common = pmap_chr(
        list(data, forecast, labs),
        ~ ..2 %>%
          filter(date == max(..1$date) - 7) %>%
          pivot_longer(contains("mid")) %>%
          filter(value == max(value, na.rm = TRUE)) %$%
          ..3[str_remove(name, "mid_")]
      )
    ) %>%
    ungroup()

  if(save) export(fits, here("outputs", paste0("manuscript_", Sys.Date(), ".rds")))

  return(fits)

}

## fit models, extract summaries, forecast
get_new_fits <- function(gis, epinow, variants = c("B.1.1.7", "B.1.351", "P.1"),
                         gen_time = 5.2, gen_time_sd = 1.9,
                         t_start = as.Date("2020-09-01"), t_end = Sys.Date() + 45,
                         min_n_var = 25, min_n_tot = 250, min_n_since = as.Date("2020-11-01"),
                         run_r_wt = TRUE,
                         save = TRUE,
                         w_diff_sq = seq(-0.3, 0.3, length = 50),
                         prop_range = c(0.01, 0.99)) {

  fits <- gis %>%
    group_by(report_country, pango_lineage) %>%
    mutate(pango_lineage = ifelse(n() > min_n_var, pango_lineage, "Other")) %>%
    ungroup() %>%
    mutate(
      y = factor(
        paste0("str", replace_na(match(pango_lineage, variants) + 1, 1)),
        levels = paste0("str", seq_len(length(variants) + 1))
      )
    ) %T>%
    {nms <<- setNames(c("Other", variants, "Joint"), c(levels(.$y), "tot"))} %>%
    group_by(report_country) %>%
    filter(
      sum(as.numeric(y) > 1) >= min_n_var,
      sum(date > min_n_since) >= min_n_tot
    ) %>%
    arrange(date) %>%
    add_times() %>%
    nest() %>%
    ungroup() %T>%
    {n_countries <<- nrow(.)} %>%
    mutate(
      ## fit models
      glm = map(data, fit_mglm2),
      ## ## predict from these models
      ## predict = map(glm, possibly(predict, NULL)),
      ## ## ## calculate rsq
      ## ## rsq = map_dbl(glm, get_glm_rsq),
      ## ## forecast proportions
      ## forecast = map2(glm, data, possibly(forecast_mglm, NULL), t_start = t_start, t_end = t_end),
      ## ## estimate r_wt
      ## r_wt = pmap(
      ##   list(glm, report_country, forecast, seq_len(n_countries)/n_countries),
      ##   possibly(est_r_wt, NULL),
      ##   epinow = epinow,
      ##   gen_time = gen_time,
      ##   gen_time_sd = gen_time_sd,
      ##   prop_range = prop_range,
      ##   type = "mglm",
      ##   run_r_wt = run_r_wt,
      ##   w_diff_sq = w_diff_sq
      ## ),
      ## ## extract summary statistics
      ## tidied = pmap(
      ##   list(
      ##     glm,
      ##     if(exists("r_wt"))
      ##       map(r_wt, possibly(
      ##         ~ filter(.x, w_diff == min(abs(w_diff)), r_mid_str1 > (1 - prop_range[2])) %$%
      ##           mean(r_mid_str1, na.rm = TRUE),
      ##         1
      ##       ))
      ##     else rep(1, n_countries),
      ##     forecast
      ##   ),
      ##   possibly(extract, NULL),
      ##   type = "mglm"
      ## ),
      ## ## ## ## ## calculate a range of r values given the generation time
      ## r_range = map2(
      ##   r_wt, tidied,
      ##   possibly(get_r_range, NULL),
      ##   gen_time = gen_time,
      ##   gen_time_sd = gen_time_sd,
      ##   run_r_wt = run_r_wt
      ## ),
      ## ## get pseudo R
      ## pseudo = map_dbl(glm, possibly(DescTools::PseudoR2, NA), "Nagelkerke"),
      ## ## get strain labels
      ## labs = rep(list(nms), length(glm)),
      ## ## get most common variant on the last day of data
      ## most_common = pmap_chr(
      ##   list(data, forecast, labs),
      ##   ~ ..2 %>%
      ##     filter(date == max(..1$date) - 7) %>%
      ##     pivot_longer(contains("mid")) %>%
      ##     filter(value == max(value, na.rm = TRUE)) %$%
      ##     ..3[str_remove(name, "mid_")]
      ## )
    ) %>%
    ungroup()

  if(save) export(fits, here("outputs", paste0("manuscript_", Sys.Date(), ".rds")))

  return(fits)

}

## fit glm model to presence/absence of b.1.1.7
fit_glm <- function(data) {
  glm(
    y ~ t,
    data = data,
    family = binomial(link = "logit")
  )
}

## fit glm model to presence/absence of b.1.1.7
fit_mglm <- function(data) {
  multinom(
    y ~ t,
    ## data = mutate(data, y = droplevels(y)),
    data = data,
    family = binomial(link = "logit"),
    model = TRUE,
    maxit = 1000
  )
}

## forecasting growth model
forecast_glm <- function(mod, data, t_start = as.Date("2020-09-01"), t_end = Sys.Date() + 45) {

  root <- min(data$date)
  dates <- seq.Date(t_start, t_end, by = 1)
  ind <- as.numeric(dates - root)

  pred <- predict(mod, newdata = tibble(t = ind), type = "link", se.fit = TRUE)

  tibble(
    date = dates,
    lower = mod$family$linkinv(pred$fit - (1.96 * pred$se.fit)),
    mid = mod$family$linkinv(pred$fit),
    upper = mod$family$linkinv(pred$fit + (1.96 * pred$se.fit))
  )

}

## forecasting multinomial regression
forecast_mglm <- function(mod, data, t_start = as.Date("2020-09-01"), t_end = Sys.Date() + 45) {

  root <- min(data$date)
  dates <- seq.Date(t_start, t_end, by = 1)
  ind <- as.numeric(dates - root)

  if(sum(table(data$y) > 0) <= 2) {

    str <- names(which(table(data$y) > 0))[2]

    sm <- summary(mod)
    coef <- sm$coefficients
    se <- sm$standard.errors

    get_probs <- function(coef, name = "mid") {
      y <- exp(coef[1] + coef[2]*ind)
      tibble(
        str1 = 1 - y/(1 + y),
        {{str}} := 1 - str1
      ) %>%
        rename_all(~ paste0(name, "_", .x))
    }

    out <- bind_cols(
      get_probs(coef, "mid"),
      get_probs(coef - 1.96*se, "lower"),
      get_probs(coef + 1.96*se, "upper")
    ) %>%
      mutate(date = dates) %>%
      select(date, everything())

  } else {

    ## use effect package to get confidence intervals
    eff <- Effect("t", mod, xlevels = list(t = ind))

    out <- data.frame(
      eff$model.matrix,
      eff$prob,
      eff$lower.prob,
      eff$upper.prob
    ) %>%
      as_tibble() %>%
      rename_all(
        ~ str_replace(.x, "L.prob.", "lower_") %>%
          str_replace("U.prob.", "upper_") %>%
          str_replace("prob.", "mid_")
      ) %>%
      mutate(date = dates) %>%
      select(date, everything(), -X.Intercept., -t)

  }

  diff <- setdiff(paste0("mid_", levels(data$y)), names(out))

  if(length(diff) > 0) {
    for(i in diff) {
      var_lower <- str_replace(i, "mid", "lower")
      var_upper <- str_replace(i, "mid", "upper")
      vals <- rep(0, nrow(out))
      out %<>% bind_cols(
        tibble(
        {{i}} := vals,
        {{var_lower}} := vals,
        {{var_upper}} := vals
        )
      )
    }
  }

  return(out)

}

## forecasting multinomial regression
forecast_mglm2 <- function(mod, data, t_start = as.Date("2020-09-01"), t_end = Sys.Date() + 45) {

  root <- min(data$date)
  dates <- seq.Date(t_start, t_end, by = 1)
  ind <- as.numeric(dates - root)

  min_dates <- data %>%
    group_by(y, .drop = FALSE) %>%
    summarise(date = min(date))

  comb <- expand.grid(dates = dates, str = levels(data$y)[-1])
  cf <- coef(mod)

  xlev <- bind_cols(
    comb,
    imap_dfc(
      setNames(min_dates$date, str_replace(min_dates$y, "str", "t")),
      function(x, y) {
        out <- as.numeric(comb$dates - x)
        out[comb$str != str_replace(y, "t", "str")] <- 0
        return(out)
      }
    )
  ) %>%
    arrange(dates) %>%
    mutate(
      p = map2_dbl(
        str, split(select(., t1:t10), seq_len(nrow(.))),
        function(str, levs) {
          mtch <- match(str, rownames(cf))
          intercept <- cf[mtch, 1]
          sums <- sum(cf[mtch, -1]*unlist(levs)[match(colnames(cf[,-1]), names(levs))]) + intercept
          return(exp(sums)/(1 + exp(sums)))
        }
      )
    )

  probs <- bind_rows(
    select(xlev, dates, str, p),
    group_by(xlev, dates) %>% summarise(str = "str1", p = 1 - sum(p, na.rm = TRUE))
  ) %>%
    arrange(dates)

  probs %>%
    group_by(dates) %>%
    summarise(p = sum(p, na.rm = TRUE)) %>%
    as.data.frame()

  probs %>%
    ggplot(aes(dates, p, fill = str)) +
    geom_line()



  if(sum(table(data$y) > 0) <= 2) {

    str <- names(which(table(data$y) > 0))[2]

    sm <- summary(mod)
    coef <- sm$coefficients
    se <- sm$standard.errors

    get_probs <- function(coef, name = "mid") {
      y <- exp(coef[1] + coef[2]*ind)
      tibble(
        str1 = 1 - y/(1 + y),
        {{str}} := 1 - str1
      ) %>%
        rename_all(~ paste0(name, "_", .x))
    }

    out <- bind_cols(
      get_probs(coef, "mid"),
      get_probs(coef - 1.96*se, "lower"),
      get_probs(coef + 1.96*se, "upper")
    ) %>%
      mutate(date = dates) %>%
      select(date, everything())

  } else {

    ## use effect package to get confidence intervals
    eff <- Effect(paste0("t", 1:8), mod, xlevels = xlev)

    out <- data.frame(
      eff$model.matrix,
      eff$prob,
      eff$lower.prob,
      eff$upper.prob
    ) %>%
      as_tibble() %>%
      rename_all(
        ~ str_replace(.x, "L.prob.", "lower_") %>%
          str_replace("U.prob.", "upper_") %>%
          str_replace("prob.", "mid_")
      ) %>%
      mutate(date = dates) %>%
      select(date, everything(), -X.Intercept., -t)

  }

  diff <- setdiff(paste0("mid_", levels(data$y)), names(out))

  if(length(diff) > 0) {
    for(i in diff) {
      var_lower <- str_replace(i, "mid", "lower")
      var_upper <- str_replace(i, "mid", "upper")
      vals <- rep(0, nrow(out))
      out %<>% bind_cols(
        tibble(
        {{i}} := vals,
        {{var_lower}} := vals,
        {{var_upper}} := vals
        )
      )
    }
  }

  return(out)

}

## show rt of different strains
vis_rt <- function(fits, country) {

  fits$r_wt[[which(fits$report_country == country)]] %>%
  pivot_longer(
    contains(c("r_wt", "r_voc", "r_tot")),
    names_to = c(".value", "var"),
    names_pattern = "(.+)_(.+)"
  ) %>%
  ggplot(aes(date, mid_r)) +
  geom_line(aes(color = var)) +#
  geom_ribbon(
    aes(ymin = lower_r, ymax = upper_r, fill = var),
    alpha = 0.5
  ) +
  scale_fill_brewer(
    name = "Variant",
    labels = c(tot = "Joint", voc = "B.1.1.7", wt = "Wild type"),
    palette = "Dark2"
  ) +
  scale_color_brewer(
    name = "Variant",
    labels = c(tot = "Joint", voc = "B.1.1.7", wt = "Wild type"),
    palette = "Dark2"
  ) +
  labs(x = NULL, y = expression(R[t])) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank())

}

## show rt of different strains
vis_rt_mglm <- function(fits, country, strains = fits$labs[[1]], base_size = 14) {

  strains <- names(fits$labs[[1]])[fits$labs[[1]] %in% strains]

  fits %<>% filter(report_country == country)

  fits$labs[[1]] %<>% modify_if(~ .x == "Joint", ~ "Observed")

  df <- fits$r_wt[[1]] %>%
    filter(w_diff == min(abs(w_diff))) %>%
    select(-contains("prop")) %>%
    pivot_longer(
      contains(c("r_l", "r_m", "r_u")),
      names_to = c(".value", "conf", "strain"),
      names_pattern = "(.+)_(.+)_(.+)"
    ) %>%
    pivot_wider(names_from = "conf", values_from = "r") %>%
    filter(strain %in% strains)

  df %>%
    ggplot(aes(date, mid)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_line(aes(color = strain, size = strain == "tot")) +
    geom_ribbon(
      data = filter(df, strain == "tot"),
      aes(ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    geom_label_repel(
      data = group_by(df, strain) %>% filter(date == max(date)),
      aes(color = strain, label = get_var_lab(fits$labs[[1]], remove = TRUE)[strain]),
      size = 6,
      segment.size = 0.5,
      segment.colour = "black",
      nudge_x = 20,
      direction = "y",
      max.overlaps = 40,
      show.legend = FALSE
    ) +
    scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 0.7), guide = FALSE) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month",
      limits = c(min(df$date), Sys.Date() + 30)
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      expand = c(0, 0)
    ) +
    ## scale_fill_manual(
    ##   name = "Variant",
    ##   labels = fits$labs[[1]],
    ##   values = c(get_var_pal(fits$labs[[1]]), tot = "black")
    ## ) +
    scale_color_manual(
      name = "Variant",
      labels = fits$labs[[1]],
      values = c(get_var_pal(fits$labs[[1]]), tot = "black"),
      guide = FALSE
    ) +
    labs(x = NULL, y = expression(R[t])) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = 'bottom'
    )

}

## show rt of different strains
vis_country_rt <- function(r_wt, labs, title = NULL, db_date = Sys.Date()) {

  strains <- names(labs)

  labs %<>% get_var_lab() %>% modify_if(~ .x == "Joint", ~ "Observed")

  df <- r_wt %>%
    filter(w_diff == min(abs(w_diff))) %>%
    select(-contains("prop")) %>%
    pivot_longer(
      contains(c("r_l", "r_m", "r_u")),
      names_to = c(".value", "conf", "strain"),
      names_pattern = "(.+)_(.+)_(.+)"
    ) %>%
    pivot_wider(names_from = "conf", values_from = "r") %>%
    filter(strain %in% strains)

  df %>%
    ggplot(aes(date, mid)) +
    geom_hline(yintercept = 1, linetype = 2) +
    geom_line(aes(color = strain, size = strain == "tot")) +
    geom_ribbon(
      data = filter(df, strain == "tot"),
      aes(ymin = lower, ymax = upper),
      alpha = 0.2
    ) +
    geom_label_repel(
      data = group_by(df, strain) %>% filter(date == max(date)),
      aes(color = strain, label = labs[strain]),
      size = 2.5,
      segment.size = 0.5,
      segment.colour = "black",
      nudge_x = 20,
      direction = "y",
      max.overlaps = 40,
      show.legend = FALSE
    ) +
    scale_size_manual(values = c("TRUE" = 1.5, "FALSE" = 0.7), guide = FALSE) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b\n%Y",
      limits = c(min(df$date), max(df$date) + 40)
    ) +
    ## scale_y_continuous(
    ##   limits = c(0.5, NA),
    ##   breaks = seq(0.6, 1.8, by = 0.2)
    ## ) +
    scale_color_manual(
      name = "Variant",
      labels = labs,
      values = c(get_var_pal(labs), tot = "black"),
      guide = FALSE
    ) +
    labs(
      x = NULL,
      y = expression(R[t]),
      title = title,
      caption = glue(
        "*Produced by WHO HQ COVID-19 analytics team | Data downloaded from ",
        "GISAID on {format(db_date, '%B %d %Y')}*"
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = 'bottom',
      plot.caption = element_markdown()
    )

}

## collect epinow data
get_epinow <- function() {

  new <- phifunc::pull_epinow() %>%
    mutate(
      across(
        report_country,
        ~ case_when(
          .x == "United States of America" ~ "USA",
          TRUE ~ .x
        )
      )
    ) %>%
    bind_rows(
      .,
      filter(., report_country == "Singapore") %>%
      mutate(report_country = "Hong Kong")
    )

  dates <- new %>%
    group_by(iso3) %>%
    summarise(min_date = min(date))

  old <- import(here("data/epinow_with_iso.csv")) %>%
    as_tibble() %>%
    left_join(dates, c(iso = "iso3")) %>%
    filter(date < min_date) %>%
    mutate(
      strat = as.character(strat),
      variable = "R"
    ) %>%
    rename(iso3 = iso) %>%
    left_join(phifunc::pull_pop_data(), "iso3") %>%
    select(-report_country, -min_date) %>%
    rename(report_country = country) %>%
    bind_rows(
      .,
      filter(., report_country == "Singapore") %>%
      mutate(report_country = "Hong Kong")
    )

  ## adjust vector in from so that it gradually becomes to by the end
  adjust <- function(from, to, range = 15) {
    if(range > length(from)) range <- length(from)
    ind <- seq(length(from) - range + 1, length(from))
    prop <- seq(0, 1, length = range)
    from[ind] <- (1 - prop)*from[ind] + prop*to
    return(from)
  }

  ## smooth Rt estimates when merging databases
  smooth <- function(old, new) {

    if(is.null(old)) return(old)
    if(is.null(new)) return(select(old, median:upper_90))

    earliest <- filter(new, date == min(date))

    map2_dfc(
      select(old, median:upper_90),
      select(earliest, median:upper_90),
      adjust
    )

  }

  ## extend old estimates if there is a gap
  fix_old <- function(old, new) {

    date_diff <- as.numeric(min(as.Date(new$date)) - (max(as.Date(old$date)) + 1))
    if(date_diff > 0 & !is.infinite(date_diff)) {
      old %<>% bind_rows(
        slice(old, rep(nrow(old), date_diff)) %>%
        mutate(date = max(old$date) + seq_len(date_diff))
      )
    }

    return(old)

  }

  ## smooth Rt estimates when joining them
  old <- left_join(
    nest(old, old = -report_country),
    nest(select(new, report_country, date, median:upper_90), new = -report_country)
  ) %>%
    mutate(
      old = map2(old, new, fix_old),
      adjust = map2(old, new, smooth),
      old = map(old, select, -(median:upper_90))
    ) %>%
    select(-new) %>%
    unnest(c(old, adjust))

  bind_rows(new, old) %>%
    arrange(iso3, date)

}

## return string after first occurence of spl
str_bef <- function(x, spl) {
  unlist(lapply(strsplit(x, spl, fixed = TRUE), '[', 1))
}

## return string before first occurence of spl
str_aft <- function(x, spl) {
  unlist(lapply(strsplit(x, spl, fixed = TRUE), '[', 2))
}

## get multionmial confidence ingeravls
get_multinom_ci <- function(x) {
  tb <- table(x)
  scimple_ci(unname(as.vector(tb)), 0.05, methods = "goodman") %>%
    transmute(
      strain = names(tb),
      lower = lower_limit,
      mid = unname(as.numeric(prop.table(tb))),
      upper = upper_limit
    )
}

## get the range of r values for a fixed growth_diff, r_wt but different w_voc
get_r_range <- function(r_wt, tidied,
                        gen_time = 5.2,
                        gen_time_sd = 1,
                        kappa = NULL,
                        run_r_wt = TRUE) {

  if(run_r_wt) {

    r_wt %>%
      select(-contains(c("prop", "tot"))) %>%
      pivot_longer(
        -contains(c("date", "w_diff", "str1")),
        names_to = c(".value", "conf", "strain"),
        names_pattern = "(.+)_(.+)_(.+)"
      ) %>%
      pivot_wider(names_from = "conf", values_from = "r") %>%
      mutate(
        lower_r_ratio = (lower - r_lower_str1)/r_lower_str1,
        mid_r_ratio = (mid - r_mid_str1)/r_mid_str1,
        upper_r_ratio = (upper - r_upper_str1)/r_upper_str1
      ) %>%
      group_by(strain, w_diff) %>%
      summarise(across(contains("r_ratio"), mean, na.rm = TRUE))

  } else return(NULL)

}

## calculate the difference in growth rates from rt and generation time
calc_growth_diff <- function(r_wt, r_voc, w_wt, w_voc, w_sd = 1.9) {
  kappa <- (w_sd/w_wt)^2
  (r_voc^kappa + (1 - r_wt^kappa)*w_voc/w_wt - 1)/(kappa*w_voc)
}

## look at the set of likely generation times and and reproduction numbers given
## the difference in growth rates
vis_ll_landscape <- function(growth_diff_est,
                             w_voc = seq(1, 10, length = 1000),
                             r_voc = seq(0.1, 10, length = 1000)) {

  expand.grid(
    w_voc = w_voc,
    r_voc = r_voc
  ) %>%
    mutate(
      growth_diff = map2_dbl(
        w_voc, r_voc,
        ~ calc_growth_diff(r_wt = 1, r_voc = .y, w_wt = 5, w_voc = .x, w_sd = 1)
      ),
      diff = abs(growth_diff - growth_diff_est)
    ) %>%
    ## arrange(diff) %>%
    ## slice(1:100) %>%
    ## ggplot(aes(r_voc, w_voc, color = diff)) +
    ## geom_point()
    ggplot(aes(r_voc, w_voc, fill = diff, z = diff)) +
    geom_raster() +
    geom_contour(bins = 50) +
    theme_minimal()

}

## look at the range of plausible r values given differences in generation time
vis_r_range <- function(fits, ylim = c(-0.25, 1)) {

  labs <- get_var_lab(fits$labs[[1]])

  fits %<>% unnest(r_range) %>%
    mutate(strain = fct_reorder(strain, mid_r_ratio, na.rm = TRUE))

  mn <- fits %>%
    group_by(w_diff, strain) %>%
    summarise(
      across(
        contains("r_ratio"),
        median, na.rm = TRUE
      )
    ) %>%
    pivot_longer(contains("r_ratio")) %>%
    mutate(strain = factor(strain, levels(fits$strain)))

  mn <- fits %>%
    group_by(w_diff, strain) %>%
    do(
      tibble(
        name = c("lower_r_ratio", "mid_r_ratio", "upper_r_ratio"),
        value = quantile(.$mid_r_ratio, c(0.25, 0.50, 0.75), na.rm = TRUE)
      )
    )

  to_percent <- function(x) paste0(ifelse(x <= 0, "", "+"), scales::percent(x, 1))

  fits %<>%
    ggplot(aes(w_diff)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_ribbon(
      aes(ymin = lower_r_ratio, ymax = upper_r_ratio, group = report_country),
      alpha = 0.15
    ) +
    ## geom_smooth(
    ##   data = mn,
    ##   aes(y = value, size = name),
    ##   color = 'black',
    ##   se = FALSE
    ## ) +
    geom_line(
      data = mn,
      aes(y = value, size = name),
    ) +
    scale_x_continuous(
      expand = c(0, 0),
      labels = to_percent,
      breaks = c(-0.2, 0, 0.2)
    ) +
    scale_y_continuous(
      labels = to_percent,
      expand = c(0, 0)
    ) +
    scale_size_manual(values = c(lower_r_ratio = 0.75, mid_r_ratio = 1.5, upper_r_ratio = 0.75)) +
    facet_wrap(~ strain, nrow = 1, labeller = labeller(strain = labs)) +
    coord_cartesian(ylim = ylim) +
    guides(size = FALSE) +
    labs(
      x = "Change in generation time",
      y = "Change in R"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 7)
    )

}

## A vectorised function to count the occurences of x in y
count_in <- function(x, y) {
  .f <- function(i, y) sum(y == i, na.rm = TRUE)
  return(sapply(x, .f, y))
}

## get variants in countries with > min_n sequences with > min_prop monthly
## proportion from the date since
get_voi <- function(gis,
                    min_n = 1000,
                    n_countries = 2,
                    props = c(0.01, 0.10),
                    since = as.Date("2020-10-01")
                    ) {

  unq <- unique(gis$pango_lineage)
  gis %>%
    group_by(report_country) %>%
    filter(n() > min_n) %>%
    group_by(report_country, after = date > since) %>%
    do(
      tibble(
        var = unq,
        n = count_in(unq, .$pango_lineage),
        prop = n/sum(n)
      )
    ) %>%
    group_by(report_country, var) %>%
    filter(sum(after) == 1) %>%
    filter(
      ifelse(sum(!after) == 0, TRUE, prop[!after] < props[1]),
      prop[after] > props[2],
      !var %in% c("", "None")
    ) %>%
    filter(after) %>%
    mutate(var = fct_reorder(var, prop, mean, na.rm = TRUE)) %$%
    levels(var)

}

## run meta analysis on var
calc_meta <- function(data, var) {
  metagen(
    TE = data[,var, drop = TRUE],
    seTE = data[,paste0(var, "_se"), drop = TRUE],
    studlab = data[,"report_country", drop = TRUE],
    comb.fixed = FALSE,
    comb.random = TRUE,
    ## method.tau = "SJ",
    hakn = TRUE,
    prediction = TRUE,
    sm = "MD"
  )
}

## random effects pooling
get_meta <- function(fits) {

  map2_dfr(
    fits$tidied,
    fits$report_country,
    possibly(~ mutate(.x, report_country = .y), NULL)
  ) %>%
    filter(
      conf == "mid",
      !is.na(growth_diff)
    ) %>%
    nest(data = -strain) %>%
    mutate(
      growth_diff = map(data, possibly(calc_meta, NULL), "growth_diff"),
      r_ratio_gamma = map(data, possibly(calc_meta, NULL), "r_ratio_gamma"),
      growth_diff_sm = map(growth_diff, summarise_meta),
      r_ratio_gamma_sm = map(r_ratio_gamma, summarise_meta)
    )

}

## random effects pooling
get_diff_meta <- function(fits, variants = head(fits$labs[[1]], -1)[-1]) {

  variants %>%
    expand.grid(., .) %>%
    as_tibble() %>%
    mutate(
      est = map2(Var1, Var2, ~ compare_str(fits, .x, .y)),
      meta = map(est, calc_meta, "mid"),
      meta_sm = map(meta, summarise_meta)
    )

}

## summarise results as mean weighted by number of sequences
get_weighted_mean <- function(fits) {

  map2(
    fits$tidied,
    fits$data,
    ~ left_join(.x, group_by(.y, strain = y) %>% summarise(n = n()), "strain")
  ) %>%
    bind_rows() %>%
    group_by(strain, conf) %>%
    summarise(
      across(
        start_se:r_ratio_gamma_se,
        ~ sum(.x[!is.na(.x)]*n[!is.na(.x)]/sum(n[!is.na(.x)]))
      )
    )

}

## get summary estimates from meta model
summarise_meta <- function(model) {
  sm <- summary(model)
  bind_cols(
    with(sm$random, tibble(mid = TE, lower, upper)),
    with(sm$predict, tibble(pred_lower = lower, pred_upper = upper))
  )
}

## pull covariants data
get_covariants <- function() {

  pull_cov_n501y() %>%
    mutate(
      report_country = modify_if(
        report_country,
        ~ .x == "United States of America",
        ~ "USA"
      )
    ) %>%
    filter(report_country != "USA" | week > as.Date("2020-11-01"))

}

## gisaid new format
get_gis_old <- function() {

  import(find_latest("gisaid_", where = here("data"))) %>%
    rename_all(epitrix::clean_labels) %>%
    mutate(
      date = as.Date(collection_date, format = "%Y-%m-%d"),
      dist = warp_distance(date, "week", every = 2),
      report_country = str_to_title(trimws(str_aft(location, "/"))),
      report_country = case_when(
        report_country %in% c("United States", "Usa") ~ "USA",
        report_country == "South Korea" ~ "Republic of Korea",
        report_country == "Russia" ~ "Russian Federation",
        report_country == "Czech Republic" ~ "Czechia",
        report_country == "Reunion" ~ "Réunion",
        TRUE ~ report_country
      ),
      pango_lineage = case_when(
        pango_lineage %in% c("B.1.427", "B.1.427") ~ "B.1.427/B.1.429",
        pango_lineage %in% c("B.1.324", "B.1.325") ~ "B.1.324/B.1.325",
        pango_lineage %in% c("B.1.526.1", "B.1.526.2", "B.1.526.3") ~ "B.1.526",
        pango_lineage %in% c("B.1.351.1", "B.1.351.2", "B.1.351.3") ~ "B.1.351",
        pango_lineage %in% c("P.1.1", "P.1.2", "P.1.3") ~ "P.1",
        TRUE ~ pango_lineage
      )
    ) %>%
    group_by(dist) %>%
    mutate(week = median(date)) %>%
    drop_na(date, pango_lineage) %>%
    filter(
      week > as.Date("2020-01-01"),
      !pango_lineage %in% c("None", "")
    ) %>%
    ungroup() %>%
    cSplit("location", " / ") %>%
    as_tibble() %>%
    mutate(across(contains("location"), as.character))

}

## gisaid new format
get_gis <- function() {

  import(find_latest("genomic_epi_", where = here("data"))) %>%
    rename_all(epitrix::clean_labels) %>%
    mutate(
      report_country = country_exposure,
      date = as.Date(date, format = "%Y-%m-%d"),
      ## import = report_country != country_exposure,
      dist = warp_distance(date, "week", every = 2),
      report_country = case_when(
        report_country %in% c("United States", "Usa") ~ "USA",
        report_country == "South Korea" ~ "Republic of Korea",
        report_country == "Russia" ~ "Russian Federation",
        report_country == "Czech Republic" ~ "Czechia",
        report_country == "Reunion" ~ "Réunion",
        TRUE ~ report_country
      ),
      pango_lineage = case_when(
        pango_lineage %in% c("B.1.427", "B.1.427") ~ "B.1.427/B.1.429",
        pango_lineage %in% c("B.1.324", "B.1.325") ~ "B.1.324/B.1.325",
        pango_lineage %in% c("B.1.526.1", "B.1.526.2", "B.1.526.3") ~ "B.1.526",
        pango_lineage %in% c("B.1.351.1", "B.1.351.2", "B.1.351.3") ~ "B.1.351",
        pango_lineage %in% c("P.1.1", "P.1.2", "P.1.3") ~ "P.1",
        grepl("AY", pango_lineage) ~ "B.1.617.2",
        TRUE ~ pango_lineage
      )
    ) %>%
    group_by(dist) %>%
    mutate(week = median(date)) %>%
    drop_na(date, pango_lineage) %>%
    filter(
      week > as.Date("2020-01-01"),
      !pango_lineage %in% c("None", ""),
    ) %>%
    ungroup() ## %>%
    ## cSplit("location", " / ") %>%
    ## as_tibble() %>%
    ## mutate(across(contains("location"), as.character))

}

## check which country has estimates for a given strain
has_strain <- function(fits, strain) {
  strain <- names(which(fits$labs[[1]] == strain))
  map(fits$tidied, ~ !is.na(.x$growth_diff[.x$strain == strain & .x$conf == "mid"])) %>%
    modify_if(~ length(.x) == 0, ~ FALSE) %>%
    unlist()
}

## wrapper for tidy function that works for nlev = 2
tidy_glm <- function(glm) {

  td <- tidy(glm)

  if(nrow(td) == 0 | !grepl("str", td$y.level[1])) {

    sm <- summary(glm)
    coef <- sm$coefficients
    se <- sm$standard.errors

    td <- tibble(
      y.level = sm$lev[-1],
      term = names(coef),
      estimate = coef,
      std.error = se,
      statistic = NA,
      p.value = NA
    )

  }

  return(td)

}

## combine ribbon and boxplot
vis_combined <- function(rib, boxp) {

  {boxp + rib & theme(legend.position = "bottom")} +
    plot_layout(
      ncol = 1,
      guides = "collect",
      heights = c(0.35, 0.65)
    ) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_text(size = 12))

}

## get variants that are in the top_n variants for at_least n countries in the last weeks

## B.1.258 - big in Eastern Europe but low increase in R
## B.1.177 - originally from Spain but also low increase
## B.1.36 - dominant in second Hong Kong wave but not elsewhere
## B.1.2 - big in US/Canda but low increase in R
get_top_strain <- function(gis, top_n = 3, at_least = 3, weeks = 4,
                           add = c("P.1", "P.2", "B.1.351", "B.1.429", "B.1.1.7",
                                   "B.1.466.2", "B.1.525", "B.1.1.519", "B.1.1.318", "B.1.324"),
                           rm = c("A", "B", "B.1", "B.1.1.1", "B.1.1", "B.1.243", "B.1.160",
                                  "B.1.221", "B.1.1.274", "B.1.1.153", "B.1.111", "B.1.1.420",
                                  "B.1.466.2", "B.1.2", "B.1.36", "B.1.177", "B.1.258")
                           ) {
  gis %>%
    mutate(dist_one = warp_distance(date, "week", every = weeks)) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one) %>%
    do(as.data.frame(prop.table(table(.$pango_lineage)))) %>%
    rename(strain = Var1, mid = Freq) %>%
    group_by(report_country) %>%
    filter(week_one == max(week_one)) %>%
    arrange(desc(mid)) %>%
    slice(seq_len(top_n)) %>%
    group_by(strain) %>%
    filter(n() > at_least) %>%
    ungroup() %>%
    mutate(strain = fct_reorder(as.character(strain), mid, mean, .desc = TRUE)) %$%
    levels(strain) %>%
    "c"(add) %>%
    setdiff(rm)
}

## compare Rt between strains within country
compare_str <- function(fits, str1, str2) {

  str1 <- names(fits$labs[[1]])[fits$labs[[1]] == str1]
  str2 <- names(fits$labs[[1]])[fits$labs[[1]] == str2]

  fits %>%
    transmute(
      report_country,
      diff = map_dbl(
        tidied,
        ~ .x$r_ratio_gamma[.x$strain == str1 & .x$conf == "mid"] -
          .x$r_ratio_gamma[.x$strain == str2 & .x$conf == "mid"]
      ),
      mid = map_dbl(
        tidied,
        ~ (1 + .x$r_ratio_gamma[.x$strain == str1 & .x$conf == "mid"])/
          (1 + .x$r_ratio_gamma[.x$strain == str2 & .x$conf == "mid"]) - 1
      ),
      ## se of difference of means is square root of sum of squares
      mid_se = map_dbl(
        tidied,
        ~ (sqrt(.x$r_ratio_gamma_se[.x$strain == str1 & .x$conf == "mid"]^2 +
                .x$r_ratio_gamma_se[.x$strain == str2 & .x$conf == "mid"]^2))/
          (1 + .x$r_ratio_gamma[.x$strain == str2 & .x$conf == "mid"])
      ),
      lower = mid - 1.96*mid_se,
      upper = mid + 1.96*mid_se
      ## lower = pmap_dbl(
      ##   list(tidied, diff, diff_se),
      ##   ~ (..2 - 1.96*..3)/(1 + ..1$r_ratio_gamma[..1$strain == str2 & ..1$conf == "upper"])
      ## ),
      ## upper = pmap_dbl(
      ##   list(tidied, diff, diff_se),
      ##   ~ (..2 + 1.96*..3)/(1 + ..1$r_ratio_gamma[..1$strain == str2 & ..1$conf == "lower"])
      ## ),
      ## lower = mid - 1.96*se,
      ## upper = mid + 1.96*se
    ) %>%
    drop_na()

}

## get proportion of strains over time
get_prop <- function(r_wt, r_diff, start, t = 0:1e5, mode = "additive") {

  growth <- calc_growth(r_wt = r_wt, r_diff = c(0, r_diff), mode = mode)
  start <- c(1 - sum(start), start)

  map2_dfc(setNames(start, paste0("str", seq_along(start))), growth, ~ .x*exp(.y*t)) %>%
    mutate(t = t) %>%
    rowwise() %>%
    mutate(total = sum(c_across(contains("str")))) %>%
    ungroup() %>%
    mutate(across(contains("str"), ~ .x/total)) %>%
    select(t, contains("str"))

}

## get the value in y closest to every element in x
get_closest_value <- function(x, y) map_dbl(x, ~ y[which.min(abs(.x - y))])

## extract r estimates from meta
get_r_estimate <- function(meta, fits, str) {

  meta %>%
    mutate(strain = fits$labs[[1]][strain]) %>%
    select(strain, r_ratio_gamma_sm) %>%
    unnest(r_ratio_gamma_sm) %>%
    filter(strain == str)

}

## get most common variants
get_most_common <- function(gis, countries = NULL,
                            n_var = 5, date_range = c(Sys.Date() - 40, Sys.Date())) {

  out <- gis %>%
    filter(
      date >= date_range[1],
      date <= date_range[2]
    )

  if(!is.null(countries)) out <- filter(out, report_country %in% countries)
  if(nrow(out) == 0 & !is.null(countries)) out <- filter(gis, report_country %in% countries)
  if(nrow(out) == 0 & is.null(countries)) out <- gis

  data.frame(table(out$pango_lineage)) %>%
    arrange(desc(Freq)) %>%
    slice(seq_len(min(c(n_var, nrow(.))))) %>%
    rename(pango_lineage = Var1, n = Freq) %>%
    mutate(pango_lineage = as.character(pango_lineage))

}

## look at locations proportion by variant
vis_locations <- function(gis, country, str = fits$labs[[1]], since = as.Date("2020-03-01"), min_n = 2) {

  gis %>%
    filter(
      report_country == country,
      date > since,
      pango_lineage %in% str
    ) %>%
    group_by(location_3) %>%
    mutate(location_3 = ifelse(n() > min_n, location_3, "Other")) %>%
    ungroup() %>%
    mutate(location_3 = fct_infreq(na_if(location_3, ""))) %>%
    ## drop_na(location_3) %>%
    ggplot(aes(y = pango_lineage, fill = location_3)) +
    geom_bar(position = "fill", color = 'black') +
    scale_fill_viridis_d() +
    labs(x = "Proportion", y = NULL, fill = "Location") +
    theme_minimal(base_size = 9) +
    guides(fill = guide_legend(ncol = 4, title.position = "top")) +
    theme(legend.position = 'bottom')

}

## visualise fits
vis_raw_location <- function(gis,
                             fits,
                             countries = unique(fits$report_country),
                             weeks = 2,
                             log = FALSE,
                             date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                             min_date = as.Date("2021-03-01"),
                             top_n = 30,
                             min_sequences = 500,
                             palette = "Dark2",
                             voi = NULL,
                             div = "location_3",
                             n_places = 20) {

  fits <- nest(gis, data = -report_country) %>%
    filter(map_dbl(data, nrow) > min_sequences)

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>% filter(report_country %in% countries)

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      slice(1:top_n)
  }

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    mutate(
      dist_one = warp_distance(date, "week", every = weeks),
      pango_lineage = modify_if(pango_lineage, ~ !.x %in% if(is.null(voi)) "Other" else voi, ~ "Other")
    ) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one, pango_lineage) %>%
    do(as.data.frame(prop.table(table(.[[div]], useNA = 'ifany')))) %>%
    rename(place = Var1, mid = Freq) %>%
    ungroup()

  levs <- df1 %>%
    group_by(report_country, place) %>%
    filter(week_one >= min_date) %>%
    summarise(mid = sum(mid)) %>%
    arrange(desc(mid)) %>%
    ungroup() %>%
    slice(1:n_places) %>%
    ungroup() %>%
    mutate(
      place = fct_reorder(as.character(place), mid, mean, .desc = TRUE)
    ) %$%
    c(levels(place), "Other")

  df1 %<>%
    mutate(
      place = as.character(place),
      place = replace(place, !place %in% levs & !is.na(place), "Other"),
      place = factor(place, levs)
    ) %>%
    group_by(report_country, week_one, pango_lineage, place) %>%
    summarise(mid = sum(mid)) %>%
    mutate(place = replace_na(place, "XXX")) %>%
    complete(report_country, week_one, place, fill = list(mid = 0)) %>%
    mutate(place = na_if(place, "XXX"))

  get_pal <- function(n, pal = "Dark2") {
    nmax <- brewer.pal.info$maxcolors[which(rownames(brewer.pal.info) == pal)]
    cols <- colorRampPalette(brewer.pal(nmax, pal))(n)
    return(cols)
  }

  pal <- if(length(unique(df1$place)) <= 8) scale_fill_brewer(palette = palette, na.value = "grey90")
         else scale_fill_manual(values = get_pal(length(unique(df1$place)), palette), na.value = "grey90")

  wrap <- if(is.null(voi)) NULL else facet_wrap(pango_lineage ~ ., ncol = 1)

  ggplot() +
    geom_area(
      data = df1,
      aes(week_one, y = mid, fill = place),
      color = 'black'
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month"
    ) +
    labs(x = NULL, y = "Proportion of sequences", fill = "State") +
    ## facet_wrap(~ report_country, ncol = 5) +
    pal +
    wrap +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom'
    )

}

## visualise fits
vis_raw_str <- function(gis,
                        fits,
                        countries = unique(fits$report_country),
                        weeks = 2,
                        log = FALSE,
                        date_range = c(as.Date("2020-09-01"), Sys.Date() + 10),
                        min_date = as.Date("2021-03-01"),
                        top_n = 30,
                        min_sequences = 500,
                        palette = "Dark2",
                        states = NULL,
                        show_strains = NULL,
                        n_str = 20) {

  fits <- nest(gis, data = -report_country) %>%
    filter(map_dbl(data, nrow) > min_sequences)

  ## set date range
  if(is.null(date_range)) date_range <- as.Date(c("2020-01-01", "2025-01-01"))

  ## remove nulls
  fits %<>% filter(report_country %in% countries)

  if(missing(countries)) {
    fits %<>%
      arrange(desc(map_dbl(data, nrow))) %>%
      slice(1:top_n)
  }

  ## extract empirical estimates
  df1 <- unnest(fits, data) %>%
    filter(date > date_range[1] & date < date_range[2]) %>%
    mutate(
      dist_one = warp_distance(date, "week", every = weeks),
      location_3 = modify_if(location_3, ~ !.x %in% if(is.null(states)) "Other" else states, ~ "Other")
    ) %>%
    group_by(dist_one) %>%
    mutate(week_one = median(date)) %>%
    group_by(report_country, week_one, location_3) %>%
    do(cbind(as.data.frame(prop.table(table(.$pango_lineage))), data.frame(n = as.vector(table(.$pango_lineage))))) %>%
    rename(str = Var1, mid = Freq) %>%
    ungroup()

  if(is.null(show_strains)) {
    levs <- df1 %>%
      group_by(report_country, str) %>%
      filter(week_one >= min_date) %>%
      summarise(mid = sum(mid)) %>%
      arrange(desc(mid)) %>%
      ungroup() %>%
      slice(1:n_str) %>%
      ungroup() %>%
      mutate(
        str = fct_reorder(as.character(str), mid, mean, .desc = TRUE)
      ) %$%
      c(levels(str), "Other")
  } else {
    levs <- c(show_strains, "Other")
  }

  df1 %<>%
    mutate(
      str = as.character(str),
      str = replace(str, !str %in% levs & !is.na(str), "Other"),
      str = factor(str, levs)
    ) %>%
    group_by(report_country, week_one, location_3, str) %>%
    summarise(mid = sum(mid), n = sum(n)) %>%
    mutate(str = replace_na(str, "XXX")) %>%
    complete(report_country, week_one, str, fill = list(mid = 0, n = 0)) %>%
    ungroup() %>%
    mutate(
      str = na_if(str, "XXX"),
      location_3 = factor(location_3, c(setdiff(unique(location_3), "Other"), "Other"))
    ) %>%
    group_by(location_3) %>%
    mutate(location_lab = glue("{location_3} (n = {sum(n)})")) %>%
    ungroup() %>%
    mutate(location_lab = fct_reorder(location_lab, as.numeric(location_3)))

  get_pal <- function(n, pal = "Dark2") {
    nmax <- brewer.pal.info$maxcolors[which(rownames(brewer.pal.info) == pal)]
    cols <- colorRampPalette(brewer.pal(nmax, pal))(n)
    return(cols)
  }

  pal <- if(length(unique(df1$str)) <= 8) scale_fill_brewer(palette = palette, na.value = "grey90")
         else scale_fill_manual(values = get_pal(length(unique(df1$str)), palette), na.value = "grey90")

  wrap <- if(is.null(states)) NULL else facet_wrap(location_lab ~ ., ncol = 2)

  ggplot() +
    geom_area(
      data = df1,
      aes(week_one, y = mid, fill = str),
      color = 'black'
    ) +
    scale_y_continuous(
      trans = ifelse(log, "log2", "identity"),
      expand = c(0, 0),
      labels = scales::percent
    ) +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month"
    ) +
    labs(x = NULL, y = "Proportion of sequences", fill = "Variant") +
    ## facet_wrap(~ report_country, ncol = 5) +
    pal +
    wrap +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      legend.position = 'bottom'
    )

}

## get the top states
get_top_states <- function(gis, n_states = 5, since = as.Date("2020-01-01")) {

  gis %>%
    filter(date >= since) %>%
    group_by(location_3) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    arrange(desc(n)) %$%
    location_3[seq_len(n_states)]

}

## show rt breakdown and modelled proportions over time
combine_rt_modelled <- function(fits, country) {

  fits %<>%
    filter(report_country == country) %>%
    mutate(forecast = map(forecast, ~ filter(.x, date %in% r_wt[[1]]$date)))

  rt <- vis_rt_mglm(fits, "India")
  modelled <- vis_modelled(fits, top_n = 1) + labs(y = expression(Contribution~to~R[t]))

  {rt + modelled &
     theme(legend.position = "bottom") &
     scale_x_date(
       expand = c(0, 0),
       date_labels = "%b\n%Y",
       date_breaks = "1 month",
       limits = c(min(unnest(fits, r_wt)$date), Sys.Date() + 30)
     )
  } +
    plot_layout(
      ncol = 1,
      guides = "collect",
      heights = c(0.6, 0.4)
    ) +
    theme(plot.tag = element_text(size = 12))

}

## show rt breakdown and modelled proportions over time
combine_casedeath_raw <- function(gis, fits, dat, vacc, stringency, country,
                                  title = NULL, levs = NULL,
                                  date_range = c(as.Date("2020-09-01"), Sys.Date())) {

  raw <- vis_raw(
    gis, country, weeks = 2,
    n_variants = 7,
    palette = "Set1",
    variants = levs,
    last_n_days = 50,
    date_range = c(as.Date("2020-09-01"), Sys.Date() - 20)
  ) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month",
      limits = c(as.Date("2020-09-01"), Sys.Date() - 20)
    )
    ## guides(fill = FALSE)

  scale_nolab <- scale_x_date(
    expand = c(0, 0),
    labels = NULL,
    name = NULL,
    date_breaks = "1 month",
    limits = c(as.Date("2020-09-01"), Sys.Date())
  )

  cases <- vis_casedeath(dat, country, "case_new", date_range = date_range) + scale_nolab
  deaths <- vis_casedeath(dat, country, "death_new", date_range = date_range) + scale_nolab

  vacc <- vis_vacc(vacc, country, date_range = c(as.Date("2020-09-01"), Sys.Date())) + scale_nolab

  stringency <- vis_stringency(stringency, country) + scale_nolab

  {cases + stringency + raw &
     theme_minimal(base_size = 11) &
     theme(
       legend.position = "bottom",
       panel.grid.minor = element_blank()
     )
  } +
    plot_layout(
      ncol = 1,
      ## guides = "collect",
      ## heights = c(0.2, 0.2, 0.2, 0.4)
      heights = c(0.25, 0.25, 0.5)
    ) +
    plot_annotation(
      title = title
    )

}

## show cases and deaths
vis_casedeath <- function(dat, country, var = "case_new",
                          date_range = c(as.Date("2020-09-01"), Sys.Date() + 10)) {

  dat %>%
    transmute(date = report_date, case_new, death_new, report_country) %>%
    filter(
      report_country == country,
      date > date_range[1] & date < date_range[2]
    ) %>%
    mutate(
      across(
        c(case_new, death_new),
        ~ slide_index_dbl(
          .x, .i = date,
          .f = mean,
          na.rm = TRUE,
          .before = lubridate::days(3),
          .after = lubridate::days(3)
        )
      )
    ) %>%
    pivot_longer(c(case_new, death_new)) %>%
    filter(name == var) %>%
    ggplot(aes(date, value)) +
    geom_line() +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month"
    ) +
    scale_y_continuous(labels = scales::comma, limits = c(0, NA)) +
    labs(
      x = NULL,
      y = ifelse(var == "case_new", "Reported cases", "Reported deaths")
    ) +
    guides(color = FALSE) +
    theme_minimal()

}

## show stringency index
vis_stringency <- function(stringency, country,
                           date_range = c(as.Date("2020-09-01"), Sys.Date() + 10)) {

  stringency %>%
    rename_all(epitrix::clean_labels) %>%
    filter(
      countryname == country,
      date > date_range[1] & date < date_range[2]
    ) %>%
    ggplot(aes(date, stringencyindex)) +
    geom_line() +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month"
    ) +
    scale_y_continuous(labels = scales::comma) +
    labs(
      x = NULL,
      y = "Stringency Index"
    ) +
    guides(color = FALSE) +
    theme_minimal()

}

## show vaccination
vis_vacc <- function(vacc, country, date_range = c(as.Date("2020-09-01"), Sys.Date() + 10)) {

  vacc %>%
    filter(
      report_country == country,
      date > date_range[1] & date < date_range[2]
    ) %>%
    mutate(
      across(
        c(people_vaccinated, people_fully_vaccinated),
        ~ slide_index_dbl(
          .x, .i = date,
          .f = mean,
          na.rm = TRUE,
          .before = lubridate::days(3),
          .after = lubridate::days(3)
        )
      )
    ) %>%
    pivot_longer(c(people_vaccinated, people_fully_vaccinated)) %>%
    mutate(name = fct_rev(name)) %>%
    ggplot(aes(date, value/population, fill = name, group = name)) +
    geom_area(position = 'identity') +
    scale_x_date(
      expand = c(0.01, 0),
      date_labels = "%b\n%Y",
      date_breaks = "1 month"
    ) +
    scale_y_continuous(labels = scales::percent) +
    scale_fill_manual(
      name = "Vaccine doses\nreceived",
      values = c("grey30", "grey60"),
      labels = c(people_vaccinated = "First dose", people_fully_vaccinated = "All doses")
    ) +
    labs(
      x = NULL,
      y = "Population vaccinated"
    ) +
    guides(color = FALSE) +
    theme_minimal()

}

## show incidence, seroprev and variants for Manaus
combine_manaus <- function(gis, fits, cases_manaus, sero_manaus,
                           date_range = c(as.Date("2020-03-01"), Sys.Date())) {

  raw <- vis_raw(
    gis, fits, "Brazil", weeks = 2,
    n_strains = 7, min_sequences = 0,
    palette = "Dark2",
    min_date = as.Date("2020-02-01"),
    date_range = date_range
  ) +
    scale_x_date(
      expand = c(0, 0),
      date_labels = "%b\n%Y",
      date_breaks = "2 month",
      limits = date_range
    )

  scale_nolab <- scale_x_date(
    expand = c(0, 0),
    labels = NULL,
    name = NULL,
    date_breaks = "2 month",
    limits = date_range
  )

  cases <- cases_manaus %>%
    mutate(
      across(
        c(casosNovos, obitosNovos),
        ~ slide_index_dbl(
          .x, .i = report_date,
          .f = mean,
          na.rm = TRUE,
          .before = lubridate::days(3),
          .after = lubridate::days(3)
        )
      )
    ) %>%
    ggplot(aes(report_date, casosNovos)) +
    geom_line() +
    scale_nolab +
    scale_y_continuous(labels = scales::comma) +
    labs(
      x = NULL,
      y = "Reported cases"
    ) +
    guides(color = FALSE) +
    theme_minimal()

  sero <- sero_manaus %>%
    select(-sens_spec_prev, -serorev_prev) %>%
    rename_at(vars(contains("serorev_prev")), str_replace, "serorev_prev", "adjusted") %>%
    rename_at(vars(contains("sens_spec_prev")), str_replace, "sens_spec_prev", "unadjusted") %>%
    pivot_longer(
      contains("adjusted"),
      names_to = c("model", ".value"),
      names_pattern = "(.+)_(.+)"
    ) %>%
    ggplot(aes(med_date, est/100, color = model)) +
    geom_point() +
    geom_line() +
    geom_errorbar(
      aes(
        ymin = lower/100,
        ymax = upper/100
      ),
      width = 0
    ) +
    scale_nolab +
    scale_y_continuous(labels = scales::percent) +
    scale_color_manual(
      values = c("#343477", "#AA8E39"),
      labels = str_to_title
    ) +
    labs(
      x = NULL,
      y = "Seroprevalence",
      color = "Estimate"
    ) +
    theme_minimal()

  {cases + sero + raw &
     theme_minimal(base_size = 10) &
     theme(
       panel.grid.minor = element_blank()
     )
  } +
    plot_layout(
      ncol = 1,
      heights = c(0.3, 0.3, 0.4)
    )

}

## shape data
shape_data <- function(data, report_country, variants, min_n_var) {

  data %>%
    group_by(pango_lineage) %>%
    mutate(pango_lineage = ifelse(n() > min_n_var, pango_lineage, "Other")) %>%
    ungroup() %>%
    mutate(
      report_country = report_country,
      y = factor(
        paste0("str", replace_na(match(pango_lineage, variants) + 1, 1)),
        levels = paste0("str", seq_len(length(variants) + 1))
      )
    ) %>%
    arrange(date) %>%
    slice(min(which(as.numeric(y) > 1)):nrow(.)) %>%
    mutate(t = as.numeric(date - min(date)))

}

## run main analysis
get_results <- function(gis,
                        n_var = 4,
                        update = TRUE,
                        run_r_wt = TRUE,
                        variants = NULL,
                        min_n_var = 0,
                        db_date = NULL,
                        days_forecast = 14,
                        last_n_days = 50,
                        gen_time = 5.2, gen_time_sd = 1.9,
                        date_range = c(as.Date("2020-09-01"), Sys.Date() + 45),
                        prop_range = c(0.01, 0.99)
                        ) {

  null_variants <- is.null(variants)

  if(update) {

    results <- gis %>%
      nest(data = -report_country) %>%
      mutate(
        variants = map(
          data,
          ~ if(null_variants)
              get_most_common(.x, n_var = n_var, date_range = c(Sys.Date() - last_n_days, Sys.Date()))$pango_lineage
          else variants
        ),
        data = pmap(
          list(data, report_country, variants),
          possibly(shape_data, NULL), min_n_var
        ),
        glm = map(data, possibly(fit_mglm, NULL)),
        ## predict from these models
        predict = map(glm, possibly(predict, NULL)),
        ## forecast proportions
        forecast = map2(glm, data, possibly(forecast_mglm, NULL), date_range[1], date_range[2]),
        ## estimate r_wt
        r_wt = pmap(
          list(glm, report_country, forecast, seq_len(nrow(.))/nrow(.)),
          possibly(est_r_wt, NULL),
          epinow = epinow,
          gen_time = gen_time,
          gen_time_sd = gen_time_sd,
          prop_range = prop_range,
          type = "mglm",
          run_r_wt = run_r_wt,
          w_diff_sq = 0
        ),
        ## extract summary statistics
        tidied = pmap(
          list(
            glm,
            if(run_r_wt)
              map(r_wt, possibly(
                ~ filter(.x, w_diff = min(abs(w_diff)), r_mid_str1 > (1 - prop_range[2])) %$%
                  mean(r_mid_str1, na.rm = TRUE),
                NULL
              ))
            else rep(1, nrow(.)),
            forecast
          ),
          possibly(extract, NULL),
          type = "mglm"
        ),
        ## ## ## ## calculate a range of r values given the generation time
        r_range = map2(
          r_wt, tidied,
          possibly(get_r_range, NULL),
          gen_time = gen_time,
          gen_time_sd = gen_time_sd,
          run_r_wt = run_r_wt
        ),
        ## get pseudo R
        pseudo = map_dbl(glm, possibly(DescTools::PseudoR2, NA), "Nagelkerke"),
        ## get strain labels
        labs = map2(
          variants, data,
          ~ setNames(c("Other", .x, "Joint"), c(levels(.y$y), "tot"))
        ),
        plot_split = pmap(
          list(
            data, labs, forecast,
            title = map2(
              report_country, data,
              ~ glue("{.x} ({format(nrow(.y), big.mark = ',', scientific=FALSE)} sequences)")
            )
          ),
          possibly(vis_country_proportions, NULL),
          weeks = 4, date_range = date_range, db_date = db_date,
          split = TRUE, days_forecast = days_forecast
        ),
        plot_combined = pmap(
          list(
            data, labs, forecast,
            title = map2(
              report_country, data,
              ~ glue("{.x} ({format(nrow(.y), big.mark = ',', scientific=FALSE)} sequences)")
            )
          ),
          possibly(vis_country_proportions, NULL),
          weeks = 3, date_range = date_range, db_date = db_date,
          split = FALSE, days_forecast = days_forecast
        ),
        plot_rt = pmap(
          list(
            r_wt = r_wt, labs = labs,
            title = map2(
              report_country, data,
              ~ glue("{.x} ({format(nrow(.y), big.mark = ',', scientific=FALSE)} sequences)")
            )
          ),
          possibly(vis_country_rt, NULL),
          db_date = db_date
        )
      ) %T>%
      export_proportions(
        here(glue("outputs/global/proportions{ifelse(null_variants,'_frequent', '_voi')}_{Sys.Date()}.rds"))
      ) %T>%
      export(
        here(glue("outputs/global/results{ifelse(null_variants, '_frequent', '_voi')}_{Sys.Date()}.rds"))
      )

    results %>%
      filter(!map_lgl(plot_split, is.null)) %$%
      walk2(
        list(plot_combined, plot_split, plot_rt),
        list(
          "Variant proportions (together) - ",
          "Variant proportions (separate) - ",
          "Disaggregated Rt - "
        ),
        ~ pwalk(
          list(.x, .y, report_country),
          ~ possibly(save_plot, NULL)(
            ..1, paste0(..2, " ", ..3, ".png"),
            width = 16.5, height = 10,
            folder = here("figures", ifelse(null_variants, "all_countries", "all_countries_voc"), ..3)
          )
        )
      )

  } else {

    results <- import(find_latest("results_frequent", where = here("outputs/global")))

  }

  return(results)

}

## export variant proportions by country by day
export_proportions <- function(fits, file) {

  fits %>%
    filter(!map_lgl(forecast, is.null)) %>%
    mutate(
      forecast = map2(
        forecast, labs,
        ~ pivot_longer(
          .x, -date,
          names_to = c("quantile", "variant"),
          names_pattern = "(.+)_(.+)"
        ) %>%
          mutate(variant = .y[variant])
      )
    ) %>%
    select(report_country, forecast) %>%
    unnest(forecast) %>%
    export(file)

}

## capitalise first letter
cap_first <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## vis HCF mask policy in detail
vis_mask_detail <- function(dat, dat2, type = c("raw", "proportion")) {

  type <- match.arg(type)

  dat2 %<>%
    pivot_longer(-c(Region, Policy)) %>%
    mutate(
      name = str_replace(cap_first(name), ", ", ",\n"),
      Policy = str_replace(Policy, " ", "\n")
    ) %>%
    left_join(select(dat, Region, `Country total`)) %>%
    ## group_by(Region) %>%
    mutate(
      label = scales::percent(value/`Country total`, 1),
      label = replace(label, value == 0, "")
    )

  if(type == "proportion") dat2 %<>% mutate(value = value/`Country total`)

  dat2 %>%
    ggplot(aes(Policy, value, fill = name)) +
    geom_bar(
      stat = 'identity',
      color = 'black'
    ) +
    geom_text(
      aes(label = label),
      size = 2.5,
      position = position_stack(vjust = 0.5)
    ) +
    facet_wrap(~ Region, nrow = 1, switch = "x") +
    scale_y_continuous(
      labels = if(type == "raw") waiver() else scales::percent,
      expand = c(0.01, 0)
    ) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x = NULL,
      y = ifelse(type == "raw", "Number of countries", "Proportion of countries"),
      fill = "Mask policy"
    ) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      strip.placement = "outside",
      legend.position = 'bottom'
    )

}

## vis HCF mask policy
vis_mask_policy <- function(dat, type = c("raw", "proportion")) {

  type <- match.arg(type)

  dat %>%
    select(-`Countries with mask guidance`, -`Country total`) %>%
    filter(Region != "TOTAL") %>%
    mutate(across(-Region, ~ as.numeric(str_remove_all(.x, " ")))) %>%
    pivot_longer(-Region) %>%
    mutate(
      name = replace(name, name %in% c("Not replied", "Unknown"), "No reply/\nUnknown") %>%
        factor(c("Risk based", "Universal", "No reply/\nUnknown"))
    ) %>%
    group_by(Region, name) %>%
    summarise(value = sum(value, na.rm = TRUE)) %>%
    ggplot(aes(Region, value, fill = name)) +
    geom_bar(
      stat = 'identity',
      color = 'black',
      position = ifelse(type == "raw", "stack", "fill")
    ) +
    scale_y_continuous(
      labels = if(type == "raw") waiver() else scales::percent,
      expand = c(0.01, 0)
    ) +
    scale_fill_manual(values = c(brewer.pal(4, "Dark2")[c(2, 3)], "grey70")) +
    labs(
      x = "WHO Region",
      y = ifelse(type == "raw", "Number of countries", "Proportion of countries"),
      fill = "Mask policy"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank())

}

## vis HCF mask policy
vis_distancing_policy <- function(dat3, type = c("raw", "proportion")) {

  type <- match.arg(type)

  dat3 %>%
    select(-`Countries with distance guidance`, -`Country total`) %>%
    filter(Region != "TOTAL") %>%
    pivot_longer(-Region) %>%
    mutate(
      name = str_replace(cap_first(name), "1m", "1.0m"),
      name = str_replace(name, "2m", "2.0m"),
    ) %>%
    ggplot(aes(Region, value, fill = name)) +
    geom_bar(
      stat = 'identity',
      color = 'black',
      position = ifelse(type == "raw", "stack", "fill")
    ) +
    scale_y_continuous(
      labels = if(type == "raw") waiver() else scales::percent,
      expand = c(0.01, 0)
    ) +
    scale_fill_manual(values = c(brewer.pal(4, "Dark2")[c(1:4)], "grey70")) +
    labs(
      x = "WHO Region",
      y = ifelse(type == "raw", "Number of countries", "Proportion of countries"),
      fill = "Distancing policy"
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank())

}

## add plus sign to percent
pminus_percent <- function(x, ...) {
  lab <- scales::percent(x, ...)
  ifelse(x > 0, paste0("+", lab), lab)
}

## pairwise difference in R between variants
vis_difference <- function(fits,
                           variants = head(fits$labs[[1]], -1)[-1],
                           base_size = 12,
                           relab = highlight_voc,
                           show_text = TRUE
                           ) {

  df <- variants %>%
    expand.grid(., .) %>%
    as_tibble() %>%
    mutate(
      est = map2(Var1, Var2, ~ compare_str(fits, .x, .y)),
      meta = map(est, calc_meta, "mid"),
      meta_sm = map(meta, summarise_meta)
    ) %>%
    filter(Var1 != Var2)

  df1 <- df %>%
    unnest(est) %>%
    group_by(Var1, Var2) %>%
    mutate(status = ifelse(lower > 0, "Higher", ifelse(upper < 0, "Lower", "Unclear")) %>%
           factor(c("Lower", "Unclear", "Higher"))) %>%
    arrange(status, mid) %>%
    mutate(n = row_number())

  df2 <- unnest(df, meta_sm) %>%
    mutate(status = ifelse(lower > 0, "Higher", ifelse(upper < 0, "Lower", "Unclear")) %>%
             factor(c("Lower", "Unclear", "Higher")))

  df3 <- tibble(Var1 = fct_inorder(levels(df1$Var1)), Var2 = fct_inorder(levels(df1$Var1))) %>%
    mutate(
      xmin = 0.72,
      xmax = 1.28,
      ymin = min(c(df1$lower)) - 0.08,
      ymax = max(c(df1$upper)) + 0.08
    )

  ggplot(df1) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(
      aes(x = 1, y = mid, group = n, color = status),
      size = 0.9,
      position = position_dodge(width = 0.5)
    ) +
    geom_errorbar(
      aes(x = 1, y = mid, group = n, color = status, ymin = lower, ymax = upper),
      position = position_dodge(width = 0.5),
      size = 0.8,
      width = 0
    ) +
    geom_errorbar(
      data = df2,
      aes(x = 1, y = mid, ymin = lower, ymax = upper),
      size = 0.9,
      width = 0
    ) +
    geom_point(
      data = df2,
      aes(x = 1, y = mid, fill = status),
      shape = 23,
      size = 2.5,
    ) +
    {if(show_text) geom_label(
      data = df2,
      aes(
        x = df3$xmin[1] + 0.13,
        y = df3$ymax[1] - 0.4,
        label = pminus_percent(mid, accuracy = 1),
        color = status
      ),
      label.size = NA,
      size = 6,
      show.legend = FALSE
    )} +
    geom_rect(
      data = df3,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      color = 'white', fill = 'white'
    ) +
    facet_grid(
      Var2 ~ Var1,
      labeller = labeller(Var1 = relab, Var2 = relab)
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(
      labels = pminus_percent,
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(min(df3$ymin), max(df3$ymax))) +
    scale_color_manual(values = c(Higher = "#D95F02", Lower ="#7570B3", Unclear = "grey")) +
    scale_fill_manual(
      values = c(Higher = "#D95F02", Lower ="#7570B3", Unclear = "grey"),
      guide = FALSE
    ) +
    labs(
      x = "Country",
      y = "Difference in R",
      color = NULL
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.background = element_rect(color = 'black'),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_blank(),
      ## axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "bottom",
      strip.text.x = element_markdown(),
      strip.text.y = element_markdown()
    )

}

## put VOCs in bold - use in conjunction with element_markdown()
highlight_voc <- function(x, voc = c("P.1", "B.1.1.7", "B.1.351", "B.1.617.2"), ...) {
  modify_if(
    get_var_lab(x, md = TRUE, ...),
    ~ grepl(paste0(voc, collapse = "|"), get_var_lab(.x)),
    ~ paste0("**", .x, "**")
  )
}

## get variant sequence numbers by country
get_country_table <- function(fits, type = 1) {

  if(type == 1) {
    fits %$%
      map2_dfr(data, report_country, ~ data.frame(report_country = .y, table(.x$y))) %>%
      transmute(
        report_country = glue("{report_country} ({prettyNum(Freq, big.mark = ',', scientific = FALSE)})"),
        var = fct_reorder(
          str_replace(get_var_lab(fits$labs[[1]][Var1]), "\n", " "),
          as.numeric(str_remove(Var1, "str"))
        ),
        Freq = Freq
      ) %>%
      arrange(desc(var)) %>%
      group_by(var) %>%
      arrange(var, desc(Freq)) %>%
      mutate(
        id = seq_len(n()),
        lab = glue("{var} [n = {prettyNum(sum(Freq), big.mark = ',', scientific = FALSE)}]")
      ) %>%
      ungroup() %>%
      filter(var != "Other") %>%
      select(-var) %>%
      pivot_wider(id_cols = id, names_from = lab, values_from = report_country) %>%
      select(-id) %>%
      mutate(across(everything(), function(x) modify_if(x, ~ grepl("\\(0\\)", .x), ~ ""))) %>%
      export(here("outputs/table.csv"))
  } else {

    fits %$%
      map2_dfr(data, report_country, ~ data.frame(report_country = .y, table(.x$y))) %>%
      transmute(
        report_country = report_country,
        var = fct_reorder(
          str_replace(get_var_lab(fits$labs[[1]][Var1]), "\n", " "),
          as.numeric(str_remove(Var1, "str"))
        ),
        Freq = Freq
      ) %>%
      filter(Freq != 0, var != "Other") %>%
      arrange(desc(var)) %>%
      group_by(var) %>%
      arrange(var, desc(Freq)) %>%
      gt::gt()

  }

}

## get confindence intervals for ratio of binomial draws
## https://stats.stackexchange.com/questions/21298/confidence-interval-around-the-ratio-of-two-proportions
get_ratio_confint <- function(x1 = 2719, n1 = 23971, x2 = 7801, n2 = 97404) {

  theta <- (x1/n1)/(x2/n2)
  se <- sqrt(1/x1 - 1/n1 + 1/x2 - 1/n2)
  c(
    lower = theta*exp(-1.96*se),
    mid = theta,
    upper = theta*exp(1.96*se)
  )

}

## rank fits by recent prevalence of given variant
rank_by_prevalence <- function(fits, variant) {

  str <- names(fits$labs[[1]])[match(variant, fits$labs[[1]])]

  fits %>%
    mutate(
      prevalence = pmap_dbl(
        list(data, forecast, labs),
        ~ ..2 %>%
          filter(date == max(..1$date) - 7) %>%
          pivot_longer(contains("mid")) %>%
          filter(grepl(str, name)) %>%
          pull(value)
      )
    ) %>%
    arrange(desc(prevalence))

}
