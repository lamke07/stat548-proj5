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

BayesOpt <- function(data_x, data_y, f_obj, par_bounds,
                     reg_model = ~1, cor_family = "PowerExponential",
                     n_starts = 20, max_iter = 50, improvement_threshold = 0.01){
  # Track EI calculations
  EI_results <- tibble()
  
  while(1){
    print(paste("Iteration", nrow(EI_results)))
    
    next_points <- BayesOpt_nextpoint(data_x = data_x, data_y = data_y, f_obj = f_obj,
                                      reg_model = reg_model, cor_family = cor_family,
                                      par_bounds, n_starts = n_starts)
    
    point_query <- next_points$next_point
    gasp_fit <- next_points$gasp_fit
    
    # # Initial and future GP fits
    # gasp_fit <- Fit(x = data_x, y = data_y,
    #                 reg_model = reg_model, cor_family = cor_family, random_error = FALSE)
    # 
    # # For multiple initial starts, find corresponding optimal points to query
    # optim_res <- get_best_EI(n_starts = n_starts, par_bounds = par_bounds, fit_obj = gasp_fit)
    # 
    # # Use optimal point query and add it to the list of EI
    # point_query <- optim_res %>%
    #   arrange(-EI) %>%
    #   slice(1)
    
    EI_results <- bind_rows(EI_results, point_query)
    
    # Check if we are expected to improve more than the threshold
    if(point_query$EI/min(gasp_fit$y) < improvement_threshold) {
      print(paste0("Expected improvement too small, stopping now: ", 100*round(point_query$EI/min(gasp_fit$y), 4), "%"))
      break
    }
    if(nrow(EI_results) >= max_iter) {
      print("Number of max iterations reached.")
      break
    }
    
    print(paste0("Expected improvement is ", 100*round(point_query$EI/min(gasp_fit$y), 4), "%... sampling more points"))
    
    # Append queried points and new function evaluation
    data_x <- bind_rows(data_x, select(point_query, -EI))
    data_y <- rbind(data_y, purrr::pmap_dbl(select(point_query, -EI), f_obj))
  }
  
  return(list(
    data = bind_cols(data_x, data_y),
    EI_results = EI_results,
    final_fit = gasp_fit
  ))
}

################################################################################
################################################################################
################################################################################

# Initial and future GP fits
gasp_fit_bayes <- GaSP2::BayesianFit(x = select(branin_data, -y), y = select(branin_data, y),
                                     reg_model = ~1, cor_family = "PowerExponential", method = "Full")
gasp_fit_bayes2 <- GaSP2::BayesianFit(x = select(branin_data, -y), y = select(branin_data, y),
                                      reg_model = ~1, cor_family = "PowerExponential", method = "Hybrid")
gasp_fit <- GaSP::Fit(x = select(branin_data, -y), y = select(branin_data, y),
                      reg_model = ~1, cor_family = "PowerExponential", random_error = FALSE)

y_pred <- GaSP::Predict(gasp_fit, x_pred, generate_coefficients = FALSE)$y_pred

par_x <- c(0, 3)
names(par_x) <- c("x1", "x2")


# Create appropriate data frame
x_pred <- data.frame(t(par_x))
colnames(x_pred) <- paste0("x", 1:length(par_x))

GaSP2::Predict(gasp_fit_bayes$GaSP_Model, data.frame(t(par_x)), generate_coefficients = TRUE, bayesian_predict = FALSE)
GaSP2::Predict(gasp_fit_bayes$GaSP_Model, x_pred, generate_coefficients = FALSE, bayesian_predict = TRUE)
GaSP2::Predict(gasp_fit_bayes$GaSP_Model, x_pred, generate_coefficients = FALSE, bayesian_predict = FALSE)
GaSP2::Predict(gasp_fit_bayes2$GaSP_Model, data.frame(t(par_x)), generate_coefficients = FALSE, bayesian_predict = FALSE)
GaSP2::Predict(gasp_fit_bayes2$GaSP_Model, data.frame(t(par_x)), generate_coefficients = FALSE, bayesian_predict = TRUE)

GaSP::Predict(gasp_fit, data.frame(t(par_x)), generate_coefficients = FALSE)
GaSP2::Predict(gasp_fit, data.frame(t(par_x)), generate_coefficients = FALSE, bayesian_predict = TRUE)
# GaSP2::Predict(gasp_fit, data.frame(t(par_x)), generate_coefficients = FALSE, bayesian_predict = FALSE)



y_pred <- GaSP2::Predict(gasp_fit$GaSP_Model, data.frame(t(par_x)), generate_coefficients = FALSE, bayesian_predict = FALSE)

?GaSP2::Predict()

names(gasp_fit)
names(gasp_fit$GaSP_Model)


