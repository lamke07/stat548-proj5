################################################################################
# Script to explore functionality of our BayesOpt algorithm (BayesOpt, BayesOpt_nextpoint)
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
par_bounds_x <- rbind(c(-5, 10),
                      c(0, 15)) %>%
  data.frame()
colnames(par_bounds_x) <- c("lower", "upper")

# Frequentist version
bayesopt_res <- BayesOpt(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
  reg_model = ~1, cor_family = "PowerExponential",
  par_bounds = par_bounds, n_starts = 100, improvement_threshold = 0.01
)

print(min(bayesopt_res$data$y)/branin(9.42478, 2.475) - 1)
print(min(bayesopt_res$data$y))

# Bayesian version
bayesopt_res_bayes <- BayesOpt(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
  reg_model = ~1, cor_family = "Matern",
  bayesian_fit = TRUE, bayesian_method = "Hybrid", bayesian_predict = FALSE,
  par_bounds = par_bounds, n_starts = 10, improvement_threshold = 0.000001
)

print(min(bayesopt_res_bayes$data$y)/branin(9.42478, 2.475) - 1)
print(min(bayesopt_res_bayes$data$y))

# Compute the next point
next_point <- BayesOpt_nextpoint(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),
  reg_model = ~1, cor_family = "PowerExponential",
  bayesian_fit = TRUE, bayesian_method = "Full", bayesian_predict = TRUE,
  par_bounds = par_bounds, n_starts = 10
)
################################################################################
################################################################################
################################################################################
# Bounds for the fried function inputs
par_bounds <- rbind(c(0,1),
                    c(0,1),
                    c(0,1),
                    c(0,1),
                    c(0,1)) %>%
  data.frame()
colnames(par_bounds) <- c("lower", "upper")

bayesopt_res <- BayesOpt(
  data_x = select(fried_data, -y), data_y = select(fried_data, y),  f_obj = fried,
  reg_model = ~1, cor_family = "PowerExponential",
  par_bounds = par_bounds, n_starts = 100, improvement_threshold = 0.000001
)

bayesopt_res$data %>%
  print(n = Inf)

################################################################################
par_bounds <- rbind(c(0,10)) %>%
  data.frame()
colnames(par_bounds) <- c("lower", "upper")

bayesopt_res <- BayesOpt(
  data_x = select(hig02_data, -y), data_y = select(hig02_data, y),  f_obj = hig02,
  reg_model = ~1, cor_family = "PowerExponential",
  par_bounds = par_bounds, n_starts = 100, improvement_threshold = 0.000001
)

bayesopt_res$data %>% 
  print(n = Inf)

################################################################################
par_bounds <- rbind(c(0,10)) %>%
  data.frame()

colnames(par_bounds) <- c("lower", "upper")

bayesopt_res <- BayesOpt(
  data_x = select(holsetal13sin_data, -y), data_y = select(holsetal13sin_data, y),  f_obj = holsetal13sin,
  reg_model = ~1, cor_family = "Matern",
  par_bounds = par_bounds, n_starts = 100, improvement_threshold = 0.000001
)

bayesopt_res_bayes <- BayesOpt(
  data_x = select(holsetal13sin_data, -y), data_y = select(holsetal13sin_data, y),  f_obj = holsetal13sin,
  reg_model = ~1, cor_family = "PowerExponential",
  bayesian_fit = TRUE, bayesian_method = "Hybrid", bayesian_predict = FALSE,
  par_bounds = par_bounds, n_starts = 10, improvement_threshold = 0.000001
)

bayesopt_res$data %>% print(n = Inf)
bayesopt_res$EI_results

################################################################################

