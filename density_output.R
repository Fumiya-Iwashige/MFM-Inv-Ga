density_estimate = function(grid_y, History, MCMC_iteration, burn_in){
  density_matrix = matrix(0.0, length(grid_y), MCMC_iteration - burn_in)
  for (i in 1:length(grid_y)) {
    for (ii in (burn_in+1):MCMC_iteration) {
      if(History[[ii]]$M_na > 0){
        S_vector = append(History[[ii]]$S_a, History[[ii]]$S_na)
        Theta_vector = append(History[[ii]]$Theta_a, History[[ii]]$Theta_na)
        Sigma_vector = append(History[[ii]]$Sigma_a, History[[ii]]$Sigma_na)
      }else{
        S_vector = History[[ii]]$S_a
        Theta_vector = History[[ii]]$Theta_a
        Sigma_vector = History[[ii]]$Sigma_a
      }
      S_vector = S_vector / sum(S_vector)
      density_matrix[i, ii-(burn_in)] = sum(S_vector * dnorm(grid_y[i], Theta_vector, sqrt(Sigma_vector)))
    }
  }
  
  density_estimate = data.frame(
    x = grid_y,
    posterior_mean = rep(0, length(grid_y)),
    quontile_0025 = rep(0, length(grid_y)),
    quontile_0975 = rep(0, length(grid_y))
  )
  
  density_estimate$posterior_mean = apply(density_matrix, 1, mean)
  for (i in 1:length(grid_y)) {
    density_estimate$quontile_0025[i] = quantile(density_matrix[i,], 0.025)
    density_estimate$quontile_0975[i] = quantile(density_matrix[i,], 0.975)
  }
  
  return(density_estimate)
}

MCMC_output_density = function(MCMC_iteration, burn_in, History){
  result_M = rep(0, MCMC_iteration-burn_in)
  result_M_na = rep(0, MCMC_iteration-burn_in)
  result_lambda = rep(0, MCMC_iteration-burn_in)
  result_U = rep(0, MCMC_iteration-burn_in)
  result_B_0 = rep(0, MCMC_iteration-burn_in)
  result_C_0 = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] = History[[i]]$M
    result_M_na[i-burn_in] = History[[i]]$M_na
    result_lambda[i-burn_in] = History[[i]]$lambda
    result_U[i-burn_in] = History[[i]]$U
    result_B_0[i-burn_in] = History[[i]]$B_0
    result_C_0[i-burn_in] = History[[i]]$C_0
  }
  
  output = c(mean(result_M),
             length(result_M_na[result_M_na==0]) / (MCMC_iteration-burn_in),
             mean(result_lambda),
             mean(result_U),
             mean(result_B_0),
             mean(result_C_0))
  
  return(output)
}

PP_output_density = function(MCMC_iteration, burn_in, History, len, S){
  result_prob = rep(0, len)
  result = rep(0, MCMC_iteration-burn_in)
  
  if(S=="M"){
    for (i in (burn_in+1):MCMC_iteration) {
      result[i-burn_in] = History[[i]]$M
    } 
  }else if(S=="k"){
    for (i in (burn_in+1):MCMC_iteration) {
      result[i-burn_in] = History[[i]]$k
    }
  }
  
  for (i in 1:len) {
    if(i < len){
      result_prob[i] = length(result[result == i]) / (MCMC_iteration-burn_in)
    }else{
      result_prob[i] = length(result[result >= len]) / (MCMC_iteration-burn_in)
    }
  }
  return(result_prob)
}