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


#Fit finite mixture model with 2 clusters for the affine model
fit_mix_model <- function(par_est_data){
  luster_1 <- par_est_data %>%
    filter(beta_0 > mean(par_est_data$beta_0))
  cluster_2 <- par_est_data %>%
    filter(beta_0 <= mean(par_est_data$beta_0))

  mu <- list( #Means from true parameters and extreme estimates
    cluster_1 = c(mean(cluster_1$beta_0),
                  mean(cluster_1$beta_1)),
    cluster_2 = c(mean(cluster_2$beta_0),
                  mean(cluster_2$beta_1))
  )

  #Fit multivariate normal finite mixture model to the estimates
  mix_model <- mvnormalmixEM(x = par_est_data[,c(3,4)], mu = mu)

  dist_table <- tibble( #Data to calculate distance
    b_0 = c(mix_model$mu[[2]][1],
            mix_model$mu[[1]][1]),
    b_1 = c(mix_model$mu[[2]][2],
            mix_model$mu[[1]][2])
  )

  return(list(
      mix_model = mix_model,
      dist = dist(dist_table)
    )
  )
}

#Plotting clusters, 1 and 2 dimensions
plot_mix_model_clusters <- function(mix_model, par_est_data){

}
