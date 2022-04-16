library(MaxPro)
library(here)
library(purrr)
library(tidyverse)

rm(list = ls())

source(here("R", "0a_functions_utils.R"))
source(here("R", "0b_functions_data.R"))


################################################################################
# LHD for the branin function
x <- MaxProLHD(n = 21, p = 2, s = 2)$Design

par <- data.frame(x1 = rescale(x[,1], target_from = -5, target_to = 10),
                  x2 = rescale(x[,2], target_from = 0, target_to = 15))
plot(par)

branin_data <- tibble(x1 = par$x1,
                      x2 = par$x2,
                      y = purrr::pmap_dbl(par, branin))

write_csv(branin_data, here("data", "branin.csv"))
################################################################################
branin(9.42478, 2.475)
branin(pi, 2.275)
branin(-pi, 12.275)