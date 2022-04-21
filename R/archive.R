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



aggregate_percs(branin_res_freq, y_min = branin_min, perc = 0.01, no_design_points = 21)

length(aggregate_percs(branin_res_freq, y_min = branin_min, perc = 0.01, no_design_points = 21))
length(aggregate_percs(branin_res_bayes, y_min = branin_min, perc = 0.01, no_design_points = 21))

purrr::map(branin_res_freq, ~purrr::map_lgl(.x, ~is_null(.x$error))) %>%
  unlist()

branin_res_freq[[7]]

# purrr::map_lgl(purrr::map(branin_res_freq[[7]], "error"), is_null)
# 
# purrr::map(branin_res_freq, ~purrr::map_dbl(.x, ~perc_threshold(.x$result$y, y_min = branin_min, perc = 0.01, no_design_points = 21))) %>%
#   unlist()
# 
# perc_threshold(branin_res_freq[[1]][[1]]$result$y, branin_min, perc = 0.01)
# 
# purrr::map(branin_res_bayes, ~purrr::map(.x, ~.x$result$y))
# purrr::map(branin_res_freq, ~purrr::map_dbl(.x, ~perc_threshold(.x$result$improvement_tracker)))
# 
# purrr::map(branin_res_freq, ~round(.x$result$improvement_tracker, 6))
# 
# perc <- 0.01
# x1 <- branin_res_freq[[1]][[1]]$result$improvement_tracker
# x <- branin_res_freq[[1]][[2]]$result$y
# # abs((min(x) - branin_min)/branin_min)
# 
# 
# 


branin_res %>%
  group_by(cor_family, bayesian_fit, bayesian_method) %>%
  mutate(no_perc_ind = sum(perc_ind > 0, na.rm = TRUE),
         no_runs = n()) %>%
  summarise(across(c(no_runs, eval_stop, eval_stop_full, min_y_stop, min_y_stop_perc, no_perc_ind, perc_ind, eval_time), median, na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(cor_family = cor_family,
            bayesian_fit = bayesian_fit,
            bayesian_method = bayesian_method,
            no_runs = paste0(no_runs_median, " (", no_runs_q0.25, " - ", no_runs_q0.75, ")"),
            eval_stop = paste0(eval_stop_median, " (", eval_stop_q0.25, " - ", eval_stop_q0.75, ")"),
            eval_stop_full = paste0(eval_stop_full_median, " (", eval_stop_full_q0.25, " - ", eval_stop_full_q0.75, ")"),
            no_perc_ind = paste0(no_perc_ind_median, " (", no_perc_ind_q0.25, " - ", no_perc_ind_q0.75, ")"),
            perc_ind = paste0(perc_ind_median, " (", perc_ind_q0.25, " - ", perc_ind_q0.75, ")"),
            eval_time = paste0(scales::comma(eval_time_median, accuracy = 0.1), " (", scales::comma(eval_time_q0.25, accuracy = 0.1), " - ", scales::comma(eval_time_q0.75, accuracy = 0.1), ")"),
            min_y_stop = paste0(scales::comma(min_y_stop_median, accuracy = 0.001), " (", scales::comma(min_y_stop_q0.25, accuracy = 0.001), " - ", scales::comma(min_y_stop_q0.75, accuracy = 0.001), ")"),
            min_y_stop_perc = paste0(scales::percent(min_y_stop_perc_median, accuracy = 0.01), " (", scales::comma(min_y_stop_perc_q0.25, accuracy = 0.01), " - ", scales::comma(min_y_stop_perc_q0.75, accuracy = 0.01), ")"),
            # (c(no_runs, eval_stop, eval_stop_full, no_perc_ind, perc_ind),
            #        ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")"))
  ) %>%
  arrange(bayesian_fit)

goldpr_res %>%
  group_by(cor_family, bayesian_fit, bayesian_method) %>%
  mutate(no_perc_ind = sum(perc_ind > 0, na.rm = TRUE),
         no_runs = n()) %>%
  summarise(across(c(no_runs, eval_stop, eval_stop_full, min_y_stop, min_y_stop_perc, no_perc_ind, perc_ind, eval_time), median, na.rm = TRUE)) %>%
  ungroup() %>%
  transmute(cor_family = cor_family,
            bayesian_fit = bayesian_fit,
            bayesian_method = bayesian_method,
            no_runs = paste0(no_runs_median, " (", no_runs_q0.25, " - ", no_runs_q0.75, ")"),
            eval_stop = paste0(eval_stop_median, " (", eval_stop_q0.25, " - ", eval_stop_q0.75, ")"),
            eval_stop_full = paste0(eval_stop_full_median, " (", eval_stop_full_q0.25, " - ", eval_stop_full_q0.75, ")"),
            no_perc_ind = paste0(no_perc_ind_median, " (", no_perc_ind_q0.25, " - ", no_perc_ind_q0.75, ")"),
            perc_ind = paste0(perc_ind_median, " (", perc_ind_q0.25, " - ", perc_ind_q0.75, ")"),
            eval_time = paste0(scales::comma(eval_time_median, accuracy = 0.1), " (", scales::comma(eval_time_q0.25, accuracy = 0.1), " - ", scales::comma(eval_time_q0.75, accuracy = 0.1), ")"),
            min_y_stop = paste0(scales::comma(min_y_stop_median, accuracy = 0.001), " (", scales::comma(min_y_stop_q0.25, accuracy = 0.001), " - ", scales::comma(min_y_stop_q0.75, accuracy = 0.001), ")"),
            min_y_stop_perc = paste0(scales::percent(min_y_stop_perc_median, accuracy = 0.01), " (", scales::comma(min_y_stop_perc_q0.25, accuracy = 0.01), " - ", scales::comma(min_y_stop_perc_q0.75, accuracy = 0.01), ")"),
            # (c(no_runs, eval_stop, eval_stop_full, no_perc_ind, perc_ind),
            #        ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")"))
  ) %>%
  arrange(bayesian_fit)



################################################################################

# Produce summary tables
branin_res %>%
  group_by(cor_family, bayesian_fit, bayesian_method) %>%
  mutate(no_perc_ind = sum(perc_ind > 0, na.rm = TRUE),
         no_runs = n()) %>%
  summarise(across(c(group_id, no_runs, eval_stop, eval_stop_full, min_y_stop, min_y_stop_perc, no_perc_ind, perc_ind, eval_time), 
                   # ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")")
                   list(median = ~median(.x, na.rm = TRUE),
                        q0.25 = ~quantile(.x, .25, na.rm = TRUE),
                        q0.75 = ~quantile(.x, .75, na.rm = TRUE)
                   )
  )) %>%
  ungroup() %>%
  transmute(group_id = group_id_median,
            cor_family = cor_family,
            bayesian_fit = bayesian_fit,
            bayesian_method = bayesian_method,
            no_runs = no_runs_median,
            eval_stop = paste0(scales::comma(eval_stop_median, accuracy = 0.1), " (", scales::comma(eval_stop_q0.25, accuracy = 0.1), "-", scales::comma(eval_stop_q0.75, accuracy = 0.1), ")"),
            eval_stop_full = paste0(scales::comma(eval_stop_full_median, accuracy = 0.1), " (", scales::comma(eval_stop_full_q0.25, accuracy = 0.1), "-", scales::comma(eval_stop_full_q0.75, accuracy = 0.1), ")"),
            no_perc_ind = paste0(no_perc_ind_median, " (", no_perc_ind_q0.25, " - ", no_perc_ind_q0.75, ")"),
            perc_ind = paste0(perc_ind_median, " (", perc_ind_q0.25, "-", perc_ind_q0.75, ")"),
            eval_time = paste0(scales::comma(eval_time_median, accuracy = 0.1), " (", scales::comma(eval_time_q0.25, accuracy = 0.1), "-", scales::comma(eval_time_q0.75, accuracy = 0.1), ")"),
            min_y_stop = paste0(scales::comma(min_y_stop_median, accuracy = 0.001), " (", scales::comma(min_y_stop_q0.25, accuracy = 0.001), "-", scales::comma(min_y_stop_q0.75, accuracy = 0.001), ")"),
            min_y_stop_perc = paste0(scales::percent(min_y_stop_perc_median, accuracy = 0.1), " (", scales::percent(min_y_stop_perc_q0.25, accuracy = 0.1), "-", scales::percent(min_y_stop_perc_q0.75, accuracy = 0.1), ")"),
            # (c(no_runs, eval_stop, eval_stop_full, no_perc_ind, perc_ind),
            #        ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")"))
  ) %>%
  arrange(bayesian_fit) %>%
  mutate(bayesian_method = ifelse(bayesian_fit == FALSE, "-", bayesian_method)) %>%
  dplyr::select(Correlation = cor_family,
                # `Bayesian Fit` = bayesian_fit,
                `Bayesian Method` = bayesian_method,
                # `Type` = group_id,
                #`#successful runs (/10)` = no_runs,
                Runtime = eval_time,
                `#Eval. until first stop` = eval_stop,
                `#Eval. until second stop` = eval_stop_full,
                `Value of f at first stop` = min_y_stop,
                `% Error at first stop` = min_y_stop_perc
  ) %>%
  readr::write_csv(here("output", "branin_summary.csv"))

goldpr_res %>%
  group_by(cor_family, bayesian_fit, bayesian_method) %>%
  mutate(no_perc_ind = sum(perc_ind > 0, na.rm = TRUE),
         no_runs = n()) %>%
  summarise(across(c(group_id, no_runs, eval_stop, eval_stop_full, min_y_stop, min_y_stop_perc, no_perc_ind, perc_ind, eval_time), 
                   # ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")")
                   list(median = ~median(.x, na.rm = TRUE),
                        q0.25 = ~quantile(.x, .25, na.rm = TRUE),
                        q0.75 = ~quantile(.x, .75, na.rm = TRUE)
                   )
  )) %>%
  ungroup() %>%
  transmute(group_id = group_id_median,
            cor_family = cor_family,
            bayesian_fit = bayesian_fit,
            bayesian_method = bayesian_method,
            no_runs = no_runs_median,
            eval_stop = paste0(scales::comma(eval_stop_median, accuracy = 0.01), " (", scales::comma(eval_stop_q0.25, accuracy = 0.01), " - ", scales::comma(eval_stop_q0.75, accuracy = 0.01), ")"),
            eval_stop_full = paste0(scales::comma(eval_stop_full_median, accuracy = 0.01), " (", scales::comma(eval_stop_full_q0.25, accuracy = 0.01), " - ", scales::comma(eval_stop_full_q0.75, accuracy = 0.01), ")"),
            no_perc_ind = paste0(no_perc_ind_median, " (", no_perc_ind_q0.25, " - ", no_perc_ind_q0.75, ")"),
            perc_ind = paste0(perc_ind_median, " (", perc_ind_q0.25, " - ", perc_ind_q0.75, ")"),
            eval_time = paste0(scales::comma(eval_time_median, accuracy = 0.1), " (", scales::comma(eval_time_q0.25, accuracy = 0.1), " - ", scales::comma(eval_time_q0.75, accuracy = 0.1), ")"),
            min_y_stop = paste0(scales::comma(min_y_stop_median, accuracy = 0.001), " (", scales::comma(min_y_stop_q0.25, accuracy = 0.001), " - ", scales::comma(min_y_stop_q0.75, accuracy = 0.001), ")"),
            min_y_stop_perc = paste0(scales::percent(min_y_stop_perc_median, accuracy = 0.01), " (", scales::percent(min_y_stop_perc_q0.25, accuracy = 0.01), " - ", scales::percent(min_y_stop_perc_q0.75, accuracy = 0.01), ")"),
            # (c(no_runs, eval_stop, eval_stop_full, no_perc_ind, perc_ind),
            #        ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")"))
  ) %>%
  arrange(bayesian_fit) %>%
  dplyr::select(#Correlation = cor_family,
    #`Bayesian Fit` = bayesian_fit,
    #`Bayesian Method` = bayesian_method,
    `Type` = group_id,
    #`#successful runs (/10)` = no_runs,
    Runtime = eval_time,
    `#Evaluations until first stopping criterion` = eval_stop,
    `#Evaluations until second stopping criterion` = eval_stop_full,
    # `y_min after first stop` = min_y_stop,
    `True error after first stop` = min_y_stop_perc
  ) %>%
  readr::write_csv(here("output", "goldpr_summary.csv"))

hart3_res %>%
  group_by(cor_family, bayesian_fit, bayesian_method) %>%
  mutate(no_perc_ind = sum(perc_ind > 0, na.rm = TRUE),
         no_runs = n()) %>%
  summarise(across(c(group_id, no_runs, eval_stop, eval_stop_full, min_y_stop, min_y_stop_perc, no_perc_ind, perc_ind, eval_time), 
                   # ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")")
                   list(median = ~median(.x, na.rm = TRUE),
                        q0.25 = ~quantile(.x, .25, na.rm = TRUE),
                        q0.75 = ~quantile(.x, .75, na.rm = TRUE)
                   )
  )) %>%
  ungroup() %>%
  transmute(group_id = group_id_median,
            cor_family = cor_family,
            bayesian_fit = bayesian_fit,
            bayesian_method = bayesian_method,
            no_runs = no_runs_median,
            eval_stop = paste0(scales::comma(eval_stop_median, accuracy = 0.01), " (", scales::comma(eval_stop_q0.25, accuracy = 0.01), " - ", scales::comma(eval_stop_q0.75, accuracy = 0.01), ")"),
            eval_stop_full = paste0(scales::comma(eval_stop_full_median, accuracy = 0.01), " (", scales::comma(eval_stop_full_q0.25, accuracy = 0.01), " - ", scales::comma(eval_stop_full_q0.75, accuracy = 0.01), ")"),
            no_perc_ind = paste0(no_perc_ind_median, " (", no_perc_ind_q0.25, " - ", no_perc_ind_q0.75, ")"),
            perc_ind = paste0(perc_ind_median, " (", perc_ind_q0.25, " - ", perc_ind_q0.75, ")"),
            eval_time = paste0(scales::comma(eval_time_median, accuracy = 0.1), " (", scales::comma(eval_time_q0.25, accuracy = 0.1), " - ", scales::comma(eval_time_q0.75, accuracy = 0.1), ")"),
            min_y_stop = paste0(scales::comma(min_y_stop_median, accuracy = 0.001), " (", scales::comma(min_y_stop_q0.25, accuracy = 0.001), " - ", scales::comma(min_y_stop_q0.75, accuracy = 0.001), ")"),
            min_y_stop_perc = paste0(scales::percent(min_y_stop_perc_median, accuracy = 0.01), " (", scales::percent(min_y_stop_perc_q0.25, accuracy = 0.01), " - ", scales::percent(min_y_stop_perc_q0.75, accuracy = 0.01), ")"),
            # (c(no_runs, eval_stop, eval_stop_full, no_perc_ind, perc_ind),
            #        ~paste0(median(.x, na.rm = TRUE), " (", quantile(.x, .25, na.rm = TRUE), " - ", quantile(.x, .75, na.rm = TRUE), ")"))
  ) %>%
  arrange(bayesian_fit) %>%
  dplyr::select(#Correlation = cor_family,
    #`Bayesian Fit` = bayesian_fit,
    #`Bayesian Method` = bayesian_method,
    `Type` = group_id,
    #`#successful runs (/10)` = no_runs,
    Runtime = eval_time,
    `#Evaluations until first stopping criterion` = eval_stop,
    `#Evaluations until second stopping criterion` = eval_stop_full,
    # `y_min after first stop` = min_y_stop,
    `True error after first stop` = min_y_stop_perc
  ) %>%
  readr::write_csv(here("output", "hart3_summary.csv"))