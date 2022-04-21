library(here)
library(tidyverse)

rm(list = ls())

source(here("R", "0b_functions_data.R"))
source(here("R", "0a_functions_utils.R"))

branin_res_freq <- readRDS(here("output", "branin_res_freq.rds"))
branin_res_bayes <- readRDS(here("output", "branin_res_bayes.rds"))

goldpr_res_freq <- readRDS(here("output", "goldpr_res_freq.rds"))
# goldpr_res_bayes <- readRDS(here("output", "goldpr_res_bayes.rds"))

hart3_res_freq <- readRDS(here("output", "hart3_res_freq.rds"))
hart3_res_bayes <- readRDS(here("output", "hart3_res_bayes.rds"))

################################################################################
# Global minima of the functions

branin_min <- branin(9.42478, 2.475)
goldpr_min <- goldpr(0,-1)
hart3_min <- hart3(0.114614, 0.555649, 0.852547)
################################################################################
################################################################################
# Functions to summarize the data

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

get_summary_table <- function(df){
  df <- df %>%
    group_by(cor_family, bayesian_fit, bayesian_method) %>%
    mutate(no_perc_ind = sum(perc_ind > 0, na.rm = TRUE),
           no_runs = n()) %>%
    summarise(across(c(group_id, no_runs, eval_stop, eval_stop_full, min_y_stop, min_y_stop_perc, no_perc_ind, perc_ind, eval_time), 
                     # ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")")
                     list(median = ~median(.x, na.rm = TRUE),
                          q0.25 = ~quantile(.x, .25, na.rm = TRUE),
                          q0.75 = ~quantile(.x, .75, na.rm = TRUE)
                     )
    )) %>%
    ungroup() %>%
    transmute(group_id = group_id_median,
              cor_family = cor_family,
              bayesian_fit = bayesian_fit,
              bayesian_method = bayesian_method,
              no_runs = no_runs_median,
              eval_stop = paste0(scales::comma(eval_stop_median, accuracy = 0.1), " (", scales::comma(eval_stop_q0.25, accuracy = 0.1), "-", scales::comma(eval_stop_q0.75, accuracy = 0.1), ")"),
              eval_stop_full = paste0(scales::comma(eval_stop_full_median, accuracy = 0.1), " (", scales::comma(eval_stop_full_q0.25, accuracy = 0.1), "-", scales::comma(eval_stop_full_q0.75, accuracy = 0.1), ")"),
              no_perc_ind = paste0(no_perc_ind_median, " (", no_perc_ind_q0.25, " - ", no_perc_ind_q0.75, ")"),
              perc_ind = paste0(perc_ind_median, " (", perc_ind_q0.25, "-", perc_ind_q0.75, ")"),
              eval_time = paste0(scales::comma(eval_time_median, accuracy = 0.1), " (", scales::comma(eval_time_q0.25, accuracy = 0.1), "-", scales::comma(eval_time_q0.75, accuracy = 0.1), ")"),
              min_y_stop = paste0(scales::comma(min_y_stop_median, accuracy = 0.001), " (", scales::comma(min_y_stop_q0.25, accuracy = 0.001), "-", scales::comma(min_y_stop_q0.75, accuracy = 0.001), ")"),
              min_y_stop_perc = paste0(scales::percent(min_y_stop_perc_median, accuracy = 0.1), " (", scales::percent(min_y_stop_perc_q0.25, accuracy = 0.1), "-", scales::percent(min_y_stop_perc_q0.75, accuracy = 0.1), ")"),
              # (c(no_runs, eval_stop, eval_stop_full, no_perc_ind, perc_ind),
              #        ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")"))
    ) %>%
    arrange(bayesian_fit) %>%
    mutate(bayesian_method = ifelse(bayesian_fit == FALSE, "-", bayesian_method)) %>%
    dplyr::select(Correlation = cor_family,
                  # `Bayesian Fit` = bayesian_fit,
                  `Bayesian Method` = bayesian_method,
                  # `Type` = group_id,
                  `#Runs` = no_runs,
                  Runtime = eval_time,
                  `#Eval. first stop` = eval_stop,
                  `#Eval. second stop` = eval_stop_full,
                  `Value of f at first stop` = min_y_stop,
                  `%Error at first stop` = min_y_stop_perc
    ) 
  
  return(df)
}
################################################################################
# Format outputs into single data frames

branin_res <- bind_rows(
  purrr::map_df(branin_res_freq, ~purrr::map_df(.x, ~.x$result$out_summary)),
  purrr::map_df(branin_res_bayes, ~purrr::map_df(.x, ~.x$result$out_summary))
  ) %>%
  mutate(min_y_stop_perc = abs((min_y_stop - branin_min)/branin_min)) %>%
  relocate(min_y_stop_perc, .after = min_y_stop) %>%
  bind_cols(
    bind_rows(
      aggregate_percs(branin_res_freq, y_min = branin_min, perc = 0.01, no_design_points = 21) %>%
        filter(error_ind == FALSE),
      aggregate_percs(branin_res_bayes, y_min = branin_min, perc = 0.01, no_design_points = 21) %>%
        filter(error_ind == FALSE)
    )
  ) %>%
  group_by(bayesian_fit, cor_family, bayesian_method) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  relocate(c(min_y_stop_perc, perc_ind), .after = min_y_stop) %>%
  relocate(group_id) %>%
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
  group_by(bayesian_fit, cor_family, bayesian_method) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  relocate(c(min_y_stop_perc, perc_ind), .after = min_y_stop) %>%
  relocate(group_id) %>%
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
  group_by(bayesian_fit, cor_family, bayesian_method) %>%
  mutate(group_id = cur_group_id()) %>%
  ungroup() %>%
  relocate(c(min_y_stop_perc, perc_ind), .after = min_y_stop) %>%
  relocate(group_id) %>%
  tibble()

################################################################################
get_summary_table(branin_res) %>%
  readr::write_csv(here("output", "branin_summary.csv"))
get_summary_table(goldpr_res) %>%
  readr::write_csv(here("output", "goldpr_summary.csv"))
get_summary_table(hart3_res) %>%
  readr::write_csv(here("output", "hart3_summary.csv"))


################################################################################
# Exploratory plots (not used)

branin_res %>%
  pivot_longer(c("eval_stop", "eval_stop_full"), names_to = "eval_type", values_to = "no_evals") %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(group_id), y = no_evals, fill = eval_type)) +
  theme_light() +
  theme(legend.position = "bottom") +
  labs(x = "Group ID", y = "Number of Evaluations"#, title = "Number of Evaluations until Stopping Criterion vs. after"
       )

branin_res %>%
  ggplot() +
  geom_boxplot(aes(x = as.factor(group_id), y = min_y_stop_perc, fill = as.factor(group_id)), show.legend = FALSE) +
  theme_light() +
  theme(legend.position = "bottom") +
  labs(x = "Group ID", y = "Relative Error"#, title = "Number of Evaluations until Stopping Criterion vs. after"
  )
  # facet_wrap(~eval_type)
  


