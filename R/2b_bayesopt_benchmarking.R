################################################################################
library(GaSP)
library(here)
library(tidyverse)

rm(list = ls())

branin_data <- read_csv(here("data", "branin.csv"))
goldpr_data <- read_csv(here("data", "goldpr.csv"))
hart3_data <- read_csv(here("data", "hart3.csv"))
# hart6_data <- read_csv(here("data", "hart6.csv"))
# fried_data <- read_csv(here("data", "fried.csv"))
# hig02_data <- read_csv(here("data", "hig02.csv"))
# holsetal13sin_data <- read_csv(here("data", "holsetal13sin.csv"))

source(here("R", "0b_functions_data.R"))
source(here("R", "0a_functions_utils.R"))
source(here("R", "0c_functions_BayesOpt.R"))

dir.create(here("output"))
################################################################################
# Set of parameters for frequentist BayesOpt
par_freq_BayesOpt <- expand.grid(
  cor_family = c("PowerExponential", "Matern")) %>%
  data.frame()

# Set of parameters for Bayesian BayesOpt
par_bayesian_BayesOpt <- expand.grid(
  cor_family = c("PowerExponential", "Matern"),
  bayesian_method = c("Full", "Hybrid"),
  bayesian_predict = c(TRUE)) %>%
  data.frame()

par_bayesian_BayesOpt_split <- asplit(par_bayesian_BayesOpt, MARGIN = 2)
par_bayesian_BayesOpt_split$bayesian_predict <- as.logical(par_bayesian_BayesOpt$bayesian_predict)

################################################################################
# Benchmark functions
n_starts <- 20
max_iter <- 500
improvement_threshold <- 10^-6

# frequentist benchmark
benchmark_freq_BayesOpt <- function(data_x, data_y, f_obj, par_bounds, par_BayesOpt, run_nr){
  
  res <- purrr::pmap(asplit(par_BayesOpt, MARGIN = 2),
                     ~purrr::safely(BayesOpt_timed)(
                       data_x = data_x, data_y = data_y, f_obj = f_obj,
                       reg_model = ~1, cor_family = ..1,
                       par_bounds = par_bounds, n_starts = n_starts, max_iter = max_iter, improvement_threshold = improvement_threshold,
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
                       par_bounds = par_bounds, n_starts = n_starts, max_iter = max_iter, improvement_threshold = improvement_threshold,
                       run_nr = run_nr
                     ))
  
  return(res)
}

################################################################################
no_of_BayesOpt <- 10
################################################################################
data_select <- branin_data
f_obj <- branin

par_bounds_x <- rbind(c(-5, 10),
                      c(0, 15)) %>%
  data.frame()
colnames(par_bounds_x) <- c("lower", "upper")

branin_res_freq <- purrr::map(1:no_of_BayesOpt,
                ~benchmark_freq_BayesOpt(data_x = select(data_select, -y), data_y = select(data_select, y),  f_obj = f_obj,
                                        par_bounds = par_bounds_x, par_BayesOpt = par_freq_BayesOpt, run_nr = .x))

saveRDS(branin_res_freq, here("output", "branin_res_freq.rds"))

branin_res_bayes <- purrr::map(1:no_of_BayesOpt,
                ~benchmark_bayesian_BayesOpt(data_x = select(data_select, -y), data_y = select(data_select, y),  f_obj = f_obj,
                                        par_bounds = par_bounds_x, par_BayesOpt = par_bayesian_BayesOpt_split, run_nr = .x))

saveRDS(branin_res_bayes, here("output", "branin_res_bayes.rds"))

################################################################################
data_select <- goldpr_data
f_obj <- goldpr

par_bounds_x <- rbind(c(-2, 2),
                      c(-2, 2)) %>%
  data.frame()
colnames(par_bounds_x) <- c("lower", "upper")

goldpr_res_freq <- purrr::map(1:no_of_BayesOpt,
                       ~benchmark_freq_BayesOpt(data_x = select(data_select, -y), data_y = select(data_select, y),  f_obj = f_obj,
                                                par_bounds = par_bounds_x, par_BayesOpt = par_freq_BayesOpt, run_nr = .x))

saveRDS(goldpr_res_freq, here("output", "goldpr_res_freq.rds"))

goldpr_res_bayes <- purrr::map(1:no_of_BayesOpt,
                        ~benchmark_bayesian_BayesOpt(data_x = select(data_select, -y), data_y = select(data_select, y),  f_obj = f_obj,
                                                     par_bounds = par_bounds_x, par_BayesOpt = par_bayesian_BayesOpt_split, run_nr = .x))

saveRDS(goldpr_res_bayes, here("output", "goldpr_res_bayes.rds"))

################################################################################
data_select <- hart3_data
f_obj <- hart3

par_bounds_x <- rbind(c(0, 1),
                      c(0, 1),
                      c(0, 1)) %>%
  data.frame()
colnames(par_bounds_x) <- c("lower", "upper")

hart3_res_freq <- purrr::map(1:no_of_BayesOpt,
                       ~benchmark_freq_BayesOpt(data_x = select(data_select, -y), data_y = select(data_select, y),  f_obj = f_obj,
                                                par_bounds = par_bounds_x, par_BayesOpt = par_freq_BayesOpt, run_nr = .x))

saveRDS(hart3_res_freq, here("output", "hart3_res_freq.rds"))

hart3_res_bayes <- purrr::map(1:no_of_BayesOpt,
                        ~benchmark_bayesian_BayesOpt(data_x = select(data_select, -y), data_y = select(data_select, y),  f_obj = f_obj,
                                                     par_bounds = par_bounds_x, par_BayesOpt = par_bayesian_BayesOpt_split, run_nr = .x))

saveRDS(hart3_res_bayes, here("output", "hart3_res_bayes.rds"))