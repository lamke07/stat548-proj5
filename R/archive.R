# p_parameters = nrow(par_bounds)


par_optim <- data.frame(
  x1 = rescale(runif(n_starts), target_from = -5, target_to = 10),
  x2 = rescale(runif(n_starts), target_from = 0, target_to = 15)) %>%
  asplit(MARGIN = 1)

# colnames(df_x) <- paste0("x", 1:length(par_x))

optim_res1 <- purrr::map(par_optim, 
                         ~optim(par = .x, 
                                fn = ac_EI_x, 
                                fit_obj = gasp_fit, mode = "negative", 
                                lower = c(-5, 0), upper = c(10, 15))
)

optim_res_df <- bind_cols(
  purrr::map_df(optim_res1, "par"),
  tibble(EI = -purrr::map_dbl(optim_res1, "value"))
)

plot(df$EI)


# optim_res <- list()
# 
# for(i in (1:n_starts)){
#   optim_res[[i]] <- optim(par_optim[i,], ac_EI_x, fit_obj = gasp_fit, lower = c(-5, 0), upper = c(10, 15))
# }


# ac_EI_x(par_x = c(4,1), fit_obj = gasp_fit)
# ac_EI_x(par_x = c(9.42478, 2.475),
#         fit_obj = gasp_fit)
# ac_EI_x(par_x = c(pi, 2.275),
#         fit_obj = gasp_fit, mode = "negative")

