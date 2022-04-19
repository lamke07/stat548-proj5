library(here)
library(tidyverse)

rm(list = ls())

source(here("R", "0b_functions_data.R"))
source(here("R", "0a_functions_utils.R"))

branin_res_freq <- readRDS(here("output", "branin_res_freq.rds"))
branin_res_bayes <- readRDS(here("output", "branin_res_bayes.rds"))

goldpr_res_freq <- readRDS(here("output", "goldpr_res_freq.rds"))
# goldpr_res_bayes <- readRDS(here("output", "goldpr_res_bayes"))

hart3_res_freq <- readRDS(here("output", "hart3_res_freq.rds"))
hart3_res_bayes <- readRDS(here("output", "hart3_res_bayes.rds"))

################################################################################
# Global minima of the functions

branin_min <- branin(9.42478, 2.475)
goldpr_min <- goldpr(0,-1)
hart3_min <- hart3(0.114614, 0.555649, 0.852547)
################################################################################
################################################################################
perc_threshold <- function(x, y_min, perc = 0.01, no_design_points){
  perc_values <- abs((x - y_min)/y_min) < perc
  
  if(any(perc_values)){
    return(which(perc_values)[1] - no_design_points)
  } else{
    return(NA)
  }
}

aggregate_percs <- function(res, y_min, perc, no_design_points){
  # Indicators of when 1% absolute error is reached
  res1 <- purrr::map(res, ~purrr::map_dbl(.x, ~perc_threshold(.x$result$y, y_min = y_min, perc = perc, no_design_points = no_design_points))) %>%
    unlist()
  
  # Indicators of whether the code returned an error or not
  res2 <- purrr::map(res, ~purrr::map_lgl(.x, ~is_null(.x$error))) %>%
    unlist()
  
  return(
    tibble(perc_ind = res1, error_ind = (res2 == FALSE))
  )
}
################################################################################

branin_res <- bind_rows(
  purrr::map_df(branin_res_freq, ~purrr::map_df(.x, ~.x$result$out_summary)),
  purrr::map_df(branin_res_bayes, ~purrr::map_df(.x, ~.x$result$out_summary))
  ) %>%
  mutate(min_y_stop_perc = abs((min_y_stop - branin_min)/branin_min)) %>%
  bind_cols(
    bind_rows(
      aggregate_percs(branin_res_freq, y_min = branin_min, perc = 0.01, no_design_points = 21) %>%
        filter(error_ind == FALSE),
      aggregate_percs(branin_res_bayes, y_min = branin_min, perc = 0.01, no_design_points = 21) %>%
        filter(error_ind == FALSE)
    )
  ) %>%
  relocate(c(min_y_stop_perc, perc_ind), .after = min_y_stop) %>%
  tibble()

goldpr_res <- bind_rows(
  purrr::map_df(goldpr_res_freq, ~purrr::map_df(.x, ~.x$result$out_summary)),
  # purrr::map_df(goldpr_res_bayes, ~purrr::map_df(.x, ~.x$result$out_summary))
  ) %>%
  mutate(min_y_stop_perc = abs((min_y_stop - goldpr_min)/goldpr_min)) %>%
  relocate(min_y_stop_perc, .after = min_y_stop) %>%
  bind_cols(
    bind_rows(
      aggregate_percs(goldpr_res_freq, y_min = goldpr_min, perc = 0.01, no_design_points = 21) %>%
        filter(error_ind == FALSE),
      # aggregate_percs(goldpr_res_bayes, y_min = hart3_min, perc = 0.01, no_design_points = 33) %>%
      #   filter(error_ind == FALSE)
    )
  ) %>%
  relocate(c(min_y_stop_perc, perc_ind), .after = min_y_stop) %>%
  tibble()

hart3_res <- bind_rows(
  purrr::map_df(hart3_res_freq, ~purrr::map_df(.x, ~.x$result$out_summary)),
  purrr::map_df(hart3_res_bayes, ~purrr::map_df(.x, ~.x$result$out_summary))
  ) %>%
  mutate(min_y_stop_perc = abs((min_y_stop - hart3_min)/hart3_min)) %>%
  bind_cols(
    bind_rows(
      aggregate_percs(hart3_res_freq, y_min = hart3_min, perc = 0.01, no_design_points = 33) %>%
        filter(error_ind == FALSE),
      aggregate_percs(hart3_res_bayes, y_min = hart3_min, perc = 0.01, no_design_points = 33) %>%
        filter(error_ind == FALSE)
    )
  ) %>%
  relocate(c(min_y_stop_perc, perc_ind), .after = min_y_stop) %>%
  tibble()

################################################################################










