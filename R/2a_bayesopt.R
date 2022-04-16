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

# set.seed(1)
bayesopt_res <- BayesOpt(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),
  f_obj = branin,
  par_bounds = par_bounds,
  n_starts = 100,
  improvement_threshold = 0.000001
)

print(min(bayesopt_res$data$y))

next_point <- BayesOpt_nextpoint(
  data_x = select(branin_data, -y), data_y = select(branin_data, y),
  f_obj = branin,
  par_bounds = par_bounds,
  n_starts = 100
)



# print(min(bayesopt_res$data$y))
# bayesopt_res$EI_results
# bayesopt_res$data %>%
#   print(n = Inf)

################################################################################
################################################################################
################################################################################
