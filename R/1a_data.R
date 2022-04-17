# setwd("~/GitHub/stat548-proj5")
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


branin_data <- tibble(x1 = par$x1,
                      x2 = par$x2,
                      y = purrr::pmap_dbl(par, branin))

write_csv(branin_data, here("data", "branin.csv"))

# Optima of the branin function
branin(9.42478, 2.475)
branin(pi, 2.275)
branin(-pi, 12.275)
################################################################################
################################################################################
# LHD for the friedman function
x <- MaxProLHD(n = 51, p = 5, s = 2)$Design

par <- data.frame(x)
colnames(par) <- paste0("x", 1:5)

fried_data <- bind_cols(
  par,
  y = purrr::pmap_dbl(par, fried)
)

write_csv(fried_data, here("data", "fried.csv"))

################################################################################
# Random sampling for the hig02 function
x <- runif(11)

par <- data.frame(rescale(x, target_from = 0, target_to = 10))
colnames(par) <- "x1"

hig02_data <- bind_cols(
  par,
  y = purrr::pmap_dbl(par, hig02)
)

write_csv(hig02_data, here("data", "hig02.csv"))

################################################################################
# Random sampling for the holsetal13sin function
x <- runif(11)

par <- data.frame(rescale(x, target_from = 0, target_to = 10))
colnames(par) <- "x1"

holsetal13sin_data <- bind_cols(
  par,
  y = purrr::pmap_dbl(par, holsetal13sin)
)

write_csv(holsetal13sin_data, here("data", "holsetal13sin.csv"))
