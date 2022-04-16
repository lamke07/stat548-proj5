################################################################################
library(GaSP)
library(here)
library(tidyverse)

rm(list = ls())

branin_data <- read_csv(here("data", "branin.csv"))

source(here("R", "0b_functions_data.R"))
source(here("R", "0a_functions_utils.R"))
source(here("R", "0c_functions_BayesOpt.R"))

################################################################################
################################################################################
# x1 in [-5, 10], x2 in [0, 15].
par_bounds <- rbind(c(-5, 10),
                    c(0, 15)) %>%
  data.frame()
colnames(par_bounds) <- c("lower", "upper")

# Parameters in form of x1 x2

# set.seed(3)
bayesopt_res <- BayesOpt(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
  reg_model = ~1, cor_family = "PowerExponential",
  par_bounds = par_bounds, n_starts = 100, improvement_threshold = 0.000001
)

print(min(bayesopt_res$data$y)/branin(9.42478, 2.475) - 1)
print(min(bayesopt_res$data$y))

# set.seed(4)
bayesopt_res_bayes <- BayesOpt(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),  f_obj = branin,
  reg_model = ~1, cor_family = "Matern",
  bayesian_fit = TRUE, bayesian_method = "Hybrid", bayesian_predict = FALSE,
  par_bounds = par_bounds, n_starts = 10, improvement_threshold = 0.000001
)

print(min(bayesopt_res_bayes$data$y)/branin(9.42478, 2.475) - 1)
print(min(bayesopt_res_bayes$data$y))


next_point <- BayesOpt_nextpoint(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),
  reg_model = ~1, cor_family = "PowerExponential",
  bayesian_fit = TRUE, bayesian_method = "Full", bayesian_predict = TRUE,
  par_bounds = par_bounds, n_starts = 10
)
################################################################################
################################################################################
################################################################################

