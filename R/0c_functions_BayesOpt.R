library(here)
library(tidyverse)
library(GaSP)
# library(GaSP2)
source(here("R", "0a_functions_utils.R"))

expected_improvement <- function(ymin, yhat, yhat_se){
  # EI formula for given ymin, yhat, yhat_se
  x_norm <- (ymin - yhat)/yhat_se
  
  return((ymin - yhat)*pnorm(x_norm) + yhat_se*dnorm(x_norm))
}

acquisition_EI <- function(par_x, fit_obj, bayesian_fit = FALSE, bayesian_predict = TRUE, return_mode = "positive"){
  # par_x is a vector of input parameters to the process
  # Mode negative for the case when it is used by optim (optim optimizes to minimum by default)
  
  # Create appropriate data frame
  x_pred <- data.frame(t(par_x))
  colnames(x_pred) <- paste0("x", 1:length(par_x))
  
  # Get parameters for expected improvement
  ymin <- min(fit_obj$y)
  
  # Predict y for either the bayesian version or not
  if(bayesian_fit == TRUE){
    y_pred <- GaSP2::Predict(fit_obj, x_pred, generate_coefficients = FALSE, bayesian_predict = bayesian_predict)
  } else if(bayesian_fit == FALSE){
    suppressWarnings(
      y_pred <- GaSP::Predict(fit_obj, x_pred, generate_coefficients = FALSE)$y_pred
    )
  }
  
  # Return positive or negative expected improvement
  if(return_mode == "positive"){
    return(expected_improvement(ymin, yhat = y_pred$Pred, yhat_se = y_pred$SE))  
  } else if(return_mode == "negative"){
    return(-expected_improvement(ymin, yhat = y_pred$Pred, yhat_se = y_pred$SE))  
  }
}

get_best_EI <- function(n_starts, par_bounds, fit_obj,
                        bayesian_fit = FALSE, bayesian_predict = TRUE){
  # Function to obtain multiple optima that the function "optim" finds in order to make the optimal EI value more robust
  
  # Sample random starts for the optim function, rescale to the right range, and write it as a list object
  suppressMessages(
    optim_starts <- purrr::map2_dfc(par_bounds$lower, par_bounds$upper, ~rescale(runif(n_starts), target_from = .x, target_to = .y))
  )
  colnames(optim_starts) <- paste0("x", 1:ncol(optim_starts))
  optim_starts_list <- asplit(optim_starts, MARGIN = 1)
  
  # Optimize EI for multiple starting points
  optim_res <- purrr::map(
    optim_starts_list, 
    ~optim(par = .x,
           fn = acquisition_EI, 
           fit_obj = fit_obj, bayesian_fit = bayesian_fit, bayesian_predict = bayesian_predict, return_mode = "negative",
           lower = par_bounds$lower, upper = par_bounds$upper)
  )
  
  # Write optimal values in a data frame
  optim_res_df <- bind_cols(
    purrr::map_df(optim_res, "par"),
    tibble(EI = -purrr::map_dbl(optim_res, "value"))
  )
  
  return(optim_res_df)
}


BayesOpt_nextpoint <- function(data_x, data_y, 
                               reg_model = ~1, cor_family = "PowerExponential",
                               bayesian_fit = FALSE, bayesian_method = "Hybrid", bayesian_predict = TRUE,
                               par_bounds, n_starts = 20){
  # Recommends next point to sample from.
  
  if(bayesian_fit == TRUE){
    gasp_fit <- GaSP2::BayesianFit(x = data_x, y = data_y,
                                   reg_model = ~1, cor_family = cor_family, method = bayesian_method)$GaSP_Model
  } else{
    # Initial and future GP fits
    gasp_fit <- GaSP::Fit(x = data_x, y = data_y,
                    reg_model = ~1, cor_family = cor_family, random_error = FALSE)
  }
  
  # For multiple initial starts, find corresponding optimal points to query
  optim_res <- get_best_EI(n_starts = n_starts, 
                           bayesian_fit = bayesian_fit, bayesian_predict = bayesian_predict, 
                           par_bounds = par_bounds, fit_obj = gasp_fit)
  
  # Use optimal point query and add it to the list of EI
  point_query <- optim_res %>%
    arrange(-EI) %>%
    slice(1)
  
  return(list(
    gasp_fit = gasp_fit,
    next_point = point_query,
    optim_EI = optim_res
  ))
}


BayesOpt <- function(data_x, data_y, f_obj, par_bounds,
                     reg_model = ~1, cor_family = "PowerExponential",
                     bayesian_fit = FALSE, bayesian_method = "Hybrid", bayesian_predict = TRUE,
                     n_starts = 20, max_iter = 50, improvement_threshold = 0.01){
  # Track EI calculations and improvements
  EI_results <- tibble()
  improvement_tracker <- vector("numeric")
  
  while(1){
    print(paste("Iteration", nrow(EI_results)))
    
    # Get a selection of the next points to use
    next_points <- BayesOpt_nextpoint(data_x = data_x, data_y = data_y,
                                      reg_model = reg_model, cor_family = cor_family,
                                      bayesian_fit = bayesian_fit, bayesian_method = bayesian_method, bayesian_predict = bayesian_predict,
                                      par_bounds, n_starts = n_starts)
    
    point_query <- next_points$next_point
    gasp_fit <- next_points$gasp_fit
    EI_results <- bind_rows(EI_results, point_query)
    
    improvement_percent <- abs(point_query$EI/min(gasp_fit$y))
    improvement_tracker <- c(improvement_tracker, improvement_percent)
    
    # Check if we are expected to improve more than the threshold
    if(improvement_percent < improvement_threshold) {
      print(paste0("Expected improvement too small, stopping now: ", 100*round(improvement_percent, 4), "%"))
      break
    }
    if(nrow(EI_results) >= max_iter) {
      print("Number of max iterations reached.")
      break
    }
    
    print(paste0("Expected improvement is ", 100*round(improvement_percent, 4), "%... testing more points"))
    
    # Append queried points and new function evaluation
    data_x <- bind_rows(data_x, select(point_query, -EI))
    data_y <- rbind(data_y, purrr::pmap_dbl(select(point_query, -EI), f_obj))
  }
  
  return(list(
    no_eval = nrow(EI_results) - 1,
    data = bind_cols(data_x, data_y),
    EI_results = EI_results,
    improvement_tracker = improvement_tracker,
    final_fit = gasp_fit,
    par_setup = list(
      reg_model = deparse1(reg_model), cor_family = cor_family,
      bayesian_fit = bayesian_fit, bayesian_method = bayesian_method, bayesian_predict = bayesian_predict,
      n_starts = n_starts, max_iter = max_iter, improvement_threshold = improvement_threshold
    )
  ))
}

BayesOpt_timed <- function(data_x, data_y, f_obj, par_bounds,
                           reg_model = ~1, cor_family = "PowerExponential",
                           bayesian_fit = FALSE, bayesian_method = "Hybrid", bayesian_predict = TRUE,
                           n_starts = 20, max_iter = 50, improvement_threshold = 0.01,
                           run_nr){
  # This is a timed version of BayesOpt for benchmarking purposes
  
  res_time <- system.time(
    res <- BayesOpt(data_x, data_y, f_obj, par_bounds,
                    reg_model, cor_family,
                    bayesian_fit, bayesian_method, bayesian_predict,
                    n_starts, max_iter, improvement_threshold)
  )
  
  no_evals_stop <- which(res$improvement_tracker < 0.01)[1] - 1
  
  out_summary <- data.frame(
    run_nr = run_nr,
    eval_stop = no_evals_stop + 1,
    eval_stop_full = length(res$improvement_tracker),
    min_y_stop = min(res$final_fit$y[1:(nrow(data_x) + no_evals_stop)]),
    eval_time = res_time[3],
    reg_model = deparse1(reg_model), cor_family = cor_family,
    bayesian_fit = bayesian_fit, bayesian_method = bayesian_method, bayesian_predict = bayesian_predict,
    n_starts = n_starts, max_iter = max_iter, improvement_threshold = improvement_threshold
    )
  
  return(list(out_summary = out_summary,
              y = res$final_fit$y,
              improvement_tracker = res$improvement_tracker))
}






