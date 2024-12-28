#Helper functions for the Here Be Dragons vignette
#Fit affine model
fit_affine_model <- function(run_no, step_size, obs_data_frame,
                             int_method,
                             prior_means = c(1,1),
                             prior_sds = c(2,2)){

  print(paste0("Run: ", run_no, " step size: ", step_size,
               " Int method: ", int_method))
  #Run the model
  suppressWarnings(
    fit <- hmde_model("affine_single_ind") |>
      hmde_assign_data(n_obs = nrow(obs_data_frame),
                       y_obs = obs_data_frame$y_obs,
                       obs_index = obs_data_frame$obs_index,
                       time = obs_data_frame$time,
                       y_bar = mean(obs_data_frame$y_obs),
                       step_size = step_size,
                       int_method = int_method,
                       prior_means = prior_means,
                       prior_sds = prior_sds)  |>  #RK4
      hmde_run(chains = 1, cores = 1, iter = 2000)
  )

  #Extract parameter estimates
  ests <- hmde_extract_estimates(model = "affine_single_ind",
                                 fit = fit,
                                 input_measurement_data = obs_data_frame)

  temp <- tibble(
    run = i,
    step_size = step_size,
    beta_0 = ests$individual_data$ind_beta_0_mean,
    beta_1 = ests$individual_data$ind_beta_1_mean
  )

  return(temp)
}


#Fit finite mixture model with 2 clusters
fit_mix_model <- function(par_est_data, par_names, true_pars){
  possible_error <- par_est_data[ #Get rows with large first parmeters
    which(par_est_data[[par_names[1]]] > mean(par_est_data[[par_names[1]]]))
    ,
  ]

  #To speed up the iterative algorithm we provide some initial conditions
  error_means <- c()
  for(i in 1:length(par_names)){
    error_means[i] <- mean(possible_error[[par_names[i]]])
  }
  mu <- list( #Means from true parameters and extreme estimates
    true = true_pars,
    error = error_means
  )

  mix_model_data <- par_est_data %>%
    dplyr::select(any_of(par_names))

  #Fit multivariate normal finite mixture model to the estimates
  mix_model <- mvnormalmixEM(x = mix_model_data, mu = mu)

  #Extract values
  post_summary_vec <- c()
  for(i in 1:2){
    par_vec_temp <- c()
    for(j in 1:length(par_names)){
      par_vec_temp[j] <- mix_model$mu[[i]][j]
    }
    post_summary_vec <- c(post_summary_vec, par_vec_temp)
  }

  #Built data frame to distance between clusters
  dist_table <- tibble()
  for(i in 1:length(par_names)){
    dist_table[[par_names[i]]] <- c(mix_model$mu[[1]][i],
                                    mix_model$mu[[2]][i])
  }

  #Build summary list
  post_summary_list <- list(
    post_summary_vec = post_summary_vec,
    error_prob = step_size_mix_models[[i]]$lambda[2],
    dist = dist(dist_table)
  )

  return(list(
      mix_model = mix_model,
      post_summary_list = post_summary_list
    )
  )
}

#Plotting clusters, 1 and 2 dimensions
plot_mix_model_clusters <- function(mix_model, par_est_data){

}
