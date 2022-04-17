################################################################################
library(GaSP)
library(here)
library(tidyverse)

rm(list = ls())

branin_data <- read_csv(here("data", "branin.csv"))
fried_data <- read_csv(here("data", "fried.csv"))
hig02_data <- read_csv(here("data", "hig02.csv"))
holsetal13sin_data <- read_csv(here("data", "holsetal13sin.csv"))

source(here("R", "0b_functions_data.R"))
source(here("R", "0a_functions_utils.R"))
source(here("R", "0c_functions_BayesOpt.R"))


################################################################################
################################################################################
# Multiple data sets
# SET MAX ITER 500
# Cases
# PowerExponential/Matern
# Fit/BayesianFit (Full/Hybrid)
# BayesianPredict (True/False)

# Improvement threshold: 0.01
# nstarts: 50

# Set of parameters for frequentist BayesOpt
par_freq_BayesOpt <- expand.grid(
  cor_family = c("PowerExponential", "Matern")) %>%
  data.frame()

# Set of parameters for Bayesian BayesOpt
par_bayesian_BayesOpt <- expand.grid(
  cor_family = c("PowerExponential", "Matern"),
  bayesian_method = c("Full", "Hybrid"),
  bayesian_predict = c(TRUE, FALSE)) %>%
  data.frame()

par_bayesian_BayesOpt_split <- asplit(par_bayesian_BayesOpt, MARGIN = 2)
par_bayesian_BayesOpt_split$bayesian_predict <- as.logical(par_bayesian_BayesOpt$bayesian_predict)

benchmark_freq_BayesOpt <- function(data_x, data_y, f_obj, par_bounds, par_BayesOpt, run_nr){
  res <- purrr::pmap(asplit(par_BayesOpt, MARGIN = 2),
                     ~purrr::safely(BayesOpt_timed)(
                       data_x = data_x, data_y = data_y, f_obj = f_obj,
                       reg_model = ~1, cor_family = ..1,
                       par_bounds = par_bounds, n_starts = 100, max_iter = 500, improvement_threshold = 0.000001,
                       run_nr = run_nr
                     ))
  
  return(res)
}

benchmark_bayesian_BayesOpt <- function(data_x, data_y, f_obj, par_bounds, par_BayesOpt, run_nr){
  res <- purrr::pmap(par_BayesOpt,
                     ~purrr::safely(BayesOpt_timed)(
                       data_x = data_x, data_y = data_y, f_obj = f_obj,
                       reg_model = ~1, cor_family = ..1,
                       bayesian_fit = TRUE, bayesian_method = ..2, bayesian_predict = ..3,
                       par_bounds = par_bounds, n_starts = 20, max_iter = 500, improvement_threshold = 0.000001,
                       run_nr = run_nr
                     ))
  
  return(res)
}

################################################################################
################################################################################
par_bounds_x <- rbind(c(-5, 10),
                      c(0, 15)) %>%
  data.frame()
colnames(par_bounds_x) <- c("lower", "upper")

x1 <- purrr::map(1:5,
                ~benchmark_freq_BayesOpt(data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
                                        par_bounds = par_bounds_x, par_BayesOpt = par_freq_BayesOpt, run_nr = .x))
x2 <- purrr::map(1:5,
                ~benchmark_bayesian_BayesOpt(data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
                                        par_bounds = par_bounds_x, par_BayesOpt = par_bayesian_BayesOpt_split, run_nr = .x))



x <- benchmark_freq_BayesOpt(data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
                             par_bounds = par_bounds_x, par_BayesOpt = par_freq_BayesOpt)


purrr::map_df(x1, ~purrr::map_df(.x, ~.x$result$out_summary))

purrr::map_df(x1, ~.x$result$out_summary)

purrr::map(x, ~.x$result$y)
purrr::map(x, ~round(.x$result$improvement_tracker, 6))