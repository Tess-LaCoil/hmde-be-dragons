#Helper functions for the Here Be Dragons vignette
#Fit affine model
fit_affine_model <- function(run_no, step_size, obs_data_frame,
                             int_method,
                             prior_pars_ind_const = c(1,2),
                             prior_pars_ind_beta_1 = c(0,2)){

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
                       prior_pars_ind_const = prior_pars_ind_const,
                       prior_pars_ind_beta_1 = prior_pars_ind_beta_1)  |>  #RK4
      hmde_run(chains = 1, cores = 1, iter = 2000)
  )

  #Extract parameter estimates
  ests <- hmde_extract_estimates(fit = fit,
                                 input_measurement_data = obs_data_frame)

  temp <- tibble(
    run = i,
    step_size = step_size,
    beta_0 = ests$individual_data$ind_beta_0_mean,
    beta_1 = ests$individual_data$ind_beta_1_mean
  )

  return(temp)
}

#Fit affine model
fit_old_affine_model <- function(run_no, step_size, obs_data_frame,
                             int_method,
                             prior_pars_ind_const = c(1,2),
                             prior_pars_ind_beta_1 = c(0,2),
                             affine_model_old){

  print(paste0("Run: ", run_no, " step size: ", step_size,
               " Int method: ", int_method))
  #Run the model

  rstan_data <- hmde_model("affine_single_ind") |>
    hmde_assign_data(n_obs = nrow(obs_data_frame),
                     y_obs = obs_data_frame$y_obs,
                     obs_index = obs_data_frame$obs_index,
                     time = obs_data_frame$time,
                     y_bar = mean(obs_data_frame$y_obs),
                     step_size = step_size,
                     int_method = int_method,
                     prior_pars_ind_const = prior_pars_ind_const,
                     prior_pars_ind_beta_1 = prior_pars_ind_beta_1)
  suppressWarnings(
    fit <- sampling(affine_model_old, rstan_data,
                    chains = 1, cores = 1, iter = 2000)
  )
  fit@model_name <- "affine_single_ind"

  #Extract parameter estimates
  ests <- hmde_extract_estimates(fit = fit,
                                 input_measurement_data = obs_data_frame)

  temp <- tibble(
    run = i,
    step_size = step_size,
    beta_0 = ests$individual_data$ind_beta_0_mean,
    beta_1 = ests$individual_data$ind_beta_1_mean
  )

  return(temp)
}

opt_init_grid <- function(rstan_data, beta_0_vals, beta_1_vals,
                          verbose = FALSE, alg = "LBFGS", affine_model){
  #Grid
  inits <- expand.grid(ind_beta_0 = beta_0_vals,
                       ind_beta_1 = beta_1_vals) %>%
    mutate(ind_const = ind_beta_0 - ind_beta_1*rstan_data$y_bar)

  post_est <- tibble()
  hessian_list <- list()

  for(i in 1:nrow(inits)){
    print(paste0("Run: ", i))

    temp <- opt_affine_model(rstan_data,
                             init_vec = inits[i,],
                             run_no = i,
                             verbose = verbose,
                             alg = alg,
                             affine_model)

    print(paste0("Runtime: ", round(temp$est_tibble_temp$runtime, digits = 4), " minutes."))

    post_est <- rbind(post_est, temp$est_tibble_temp)
    hessian_list[[i]] <- temp$hessian
  }

  return(
    list(post_est = post_est,
         hessian_list = hessian_list)
  )
}

opt_affine_model <- function(rstan_data, init_vec = NULL, run_no,
                             verbose = FALSE, alg = "LBFGS", affine_model){
  if(!is.null(init_vec)){
    init <- list(ind_y_0 = rstan_data$y_obs[1],
       ind_const = init_vec$ind_const,
       ind_beta_1 = init_vec$ind_beta_1,
       y_hat = rstan_data$y_obs,
       ind_beta_0 = init_vec$ind_beta_0,
       pars = c(init_vec$ind_beta_0,
                init_vec$ind_beta_1),
       y_temp = 1)

    start_time <- Sys.time()
    fit <- optimizing(affine_model, rstan_data,
                      alg = alg,
                      hessian = TRUE,
                      init = init,
                      verbose = verbose)
    end_time <- Sys.time()

    est_tibble_temp <- tibble(
      run = run_no,
      return = fit$return_code,
      runtime = difftime(end_time, start_time, units = "mins"),
      beta_0_init = init_vec$ind_beta_0,
      beta_1_init = init_vec$ind_beta_1,
      beta_0 = fit$par[["ind_beta_0"]],
      beta_1 = fit$par[["ind_beta_1"]]
    )
  } else {

    start_time <- Sys.time()
    fit <- optimizing(affine_model, rstan_data,
                      alg = alg,
                      hessian = TRUE,
                      verbose = verbose)
    end_time <- Sys.time()

    est_tibble_temp <- tibble(
      run = run_no,
      return = fit$return_code,
      lp_val = fit$value,
      runtime = difftime(end_time, start_time, units = "mins"),
      beta_0 = fit$par[["ind_beta_0"]],
      beta_1 = fit$par[["ind_beta_1"]]
    )
  }

  return_data <- list(
    est_tibble_temp = est_tibble_temp,
    hessian = fit$hessian
  )

  return(return_data)
}

#Fit finite mixture model with 2 clusters for the affine model
fit_mix_model <- function(est_data){ #Pass in only the columns with parameter values
  cluster_1 <- est_data %>%
    filter(beta_0 <= mean(est_data$beta_0))
  cluster_2 <- est_data %>%
    filter(beta_0 > mean(est_data$beta_0))

  mu <- list( #Means from true parameters and extreme estimates
    cluster_1 = c(mean(cluster_1$beta_0),
                  mean(cluster_1$beta_1)),
    cluster_2 = c(mean(cluster_2$beta_0),
                  mean(cluster_2$beta_1))
  )

  #Fit multivariate normal finite mixture model to the estimates
  mix_model <- mvnormalmixEM(x = est_data, mu = mu)

  return(mix_model)
}

#Distance calculation
calculate_cluster_distances <- function(mix_model, centre){
  dist_table <- tibble( #Data to calculate Euclidean distance
    b_0 = c(mix_model$mu[[2]][1],
            mix_model$mu[[1]][1]),
    b_1 = c(mix_model$mu[[2]][2],
            mix_model$mu[[1]][2])
  )

  dist_vec <- c(m_dist_c1 =  mahalanobis(mix_model$mu[[1]],
              center = centre,
              cov = mix_model$sigma[[1]]),
              m_dist_c2 = mahalanobis(mix_model$mu[[2]],
                                      center = centre,
                                      cov = mix_model$sigma[[2]]),
              cluster_euclid_dist = dist(dist_table)
  )
  return(dist_vec)
}

#Plotting clusters, 1 and 2 dimensions
plot_mix_model_clusters <- function(mix_model, par_est_data){

}
