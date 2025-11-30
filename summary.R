getDahl <- function(MFMfit, burn)
{
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}

MCMC_output <- function(MCMC_iteration, burn_in, History, true_allocation, Method){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  ##        true_allocation : true labels 
  ##        Method : Method name (character), ex: "MFM-IGau"
  
  ################################################################
  
  ## Return: A data frame containing the following:
  ##        Method : Method name (character)
  ##        post.mean.ARI : posterior mean of adjusted Rand index
  ##        post.map.ARI : MAP of adjusted Rand index
  ##        post.mean.cluster : posterior mean of number of clusters k
  ##        post.medi.cluster : posterior mode or medisn of number of clusters k
  ##        post.mean.comp : posterior mean of number of components M
  ##        post.prop.M_na : posterior probability of P(M_na | data) = 0
  
  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result.ARI = rep(0, MCMC_iteration-burn_in)
  result.cluster = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result.ARI[i-burn_in] <- adj.rand.index(History[[i]]$c, true_allocation)
    result.cluster[i-burn_in] <- History[[i]]$k
  }
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  
  df_results  <- data.frame(
    Method = NA,
    post.mean.ARI = NA,
    post.map.ARI = NA,
    post.mean.cluster = NA,
    post.medi.cluster = NA,
    post.mean.comp = NA,
    post.medi.comp = NA,
    post.prop.M_na = NA
  )
  output <- data.frame(
              Method = Method,
              post.mean.ARI = mean(result.ARI),
              post.map.ARI = adj.rand.index(result_c, true_allocation),
              post.mean.cluster = mean(result.cluster),
              post.medi.cluster = median(result.cluster),
              post.mean.comp = mean(result_M),
              post.medi.comp = median(result_M),
              post.prop.M_na = length(result_M_na[result_M_na==0]) / ( MCMC_iteration-(burn_in+1) )
  )
  
  return(output)
}

MCMC_output2 <- function(MCMC_iteration, burn_in, History, true_allocation, Method){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  ##        true_allocation : true labels 
  ##        Method : Method name (character), ex: "MFM-IGau"
  
  ################################################################
  
  ## Return: A vector consisting of the following elements:
  ##        Method name (character)
  ##        posterior mean of adjusted Rand index
  ##        MAP of adjusted Rand index
  ##        posterior mean of number of clusters k
  ##        posterior mode or medisn of number of clusters k
  ##        posterior mean of number of components M
  ##        posterior probability of P(M_na | data) = 0
  ##        the size of the data used
  
  data_size <- length(true_allocation)
  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result.ARI = rep(0, MCMC_iteration-burn_in)
  result.cluster = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result.ARI[i-burn_in] <- adj.rand.index(History[[i]]$c, true_allocation)
    result.cluster[i-burn_in] <- History[[i]]$k
  }
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  
  output <- c(Method, 
              mean(result.ARI),
              adj.rand.index(result_c, true_allocation),
              mean(result.cluster),
              median(result.cluster),
              mean(result_M),
              median(result_M),
              length(result_M_na[result_M_na==0]) / ( MCMC_iteration-burn_in ),
              data_size
  )
  
  return(output)
}

MCMC_output_EP <- function(MCMC_iteration, burn_in, History, true_allocation, Method, shape){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  ##        true_allocation : true labels 
  ##        Method : Method name (character), ex: "MFM-IGau"
  ##        shape : shape parameter of h (Gamma(gamma, 1) or IGau(Alpha, 1))
  
  ################################################################
  
  ## Return: A vector consisting of the following elements:
  ##        posterior mean of adjusted Rand index
  ##        MAP of adjusted Rand index
  ##        posterior mean of number of clusters k
  ##        posterior mode or medisn of number of clusters k
  ##        posterior mean of number of components M
  ##        posterior mode or medisn of number of components M
  ##        posterior 5% quantile of M  
  ##        posterior 95% quantile of M  
  ##        posterior probability of P(M_na | data) = 0
  ##        shape parameter of h
  
  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result.ARI = rep(0, MCMC_iteration-burn_in)
  result.cluster = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result.ARI[i-burn_in] <- adj.rand.index(History[[i]]$c, true_allocation)
    result.cluster[i-burn_in] <- History[[i]]$k
  }
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  
  output <- c(mean(result.ARI),
              adj.rand.index(result_c, true_allocation),
              mean(result.cluster),
              median(result.cluster),
              mean(result_M),
              median(result_M),
              as.numeric(quantile(result_M, probs = c(0.05))),
              as.numeric(quantile(result_M, probs = c(0.95))),
              length(result_M_na[result_M_na==0]) / ( MCMC_iteration-(burn_in+1) ),
              shape
  )
  
  return(output)
}

MCMC_output_thyroid <- function(MCMC_iteration, burn_in, History, true_allocation, Method, qM){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  ##        true_allocation : true labels 
  ##        Method : Method name (character), ex: "MFM-IGau"
  ##        qM : prior distribution of M-1 (character)
  
  ################################################################
  
  ## Return: A vector consisting of the following elements:
  ##        Method name (character)
  ##        qM (character)
  ##        posterior mean of adjusted Rand index
  ##        MAP of adjusted Rand index
  ##        posterior mean of number of clusters k
  ##        posterior mode or medisn of number of clusters k
  ##        posterior 25% quantile of k  
  ##        posterior 75% quantile of k  
  ##        posterior mean of number of components M
  ##        posterior mode or medisn of number of components M
  ##        posterior 25% quantile of M  
  ##        posterior 75% quantile of M  
  ##        posterior probability of P(M_na | data) = 0

  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result.ARI = rep(0, MCMC_iteration-burn_in)
  result.cluster = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result.ARI[i-burn_in] <- adj.rand.index(History[[i]]$c, true_allocation)
    result.cluster[i-burn_in] <- History[[i]]$k
  }
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  
  output <- c(Method, 
              qM,
              mean(result.ARI),
              adj.rand.index(result_c, true_allocation),
              mean(result.cluster),
              median(result.cluster),
              as.numeric(quantile(result.cluster, probs = c(0.25))),
              as.numeric(quantile(result.cluster, probs = c(0.75))),
              mean(result_M),
              median(result_M),
              as.numeric(quantile(result_M, probs = c(0.25))),
              as.numeric(quantile(result_M, probs = c(0.75))),
              length(result_M_na[result_M_na==0]) / ( MCMC_iteration-burn_in )
  )
  
  return(output)
}

MCMC_output_commu_sim <- function(MCMC_iteration, burn_in, History, true_allocation, Ad_matrix){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  ##        true_allocation : true labels 
  ##        Ad_matrix : An adjacent matrix of size (node size × node size).

  ################################################################
  
  ## Return: A vector consisting of the following elements:
  ##        MAP of modularity
  ##        posterior mean of adjusted Rand index
  ##        MAP of adjusted Rand index
  ##        posterior mean of number of clusters k
  ##        posterior mode or medisn of number of clusters k
  ##        posterior 25% quantile of k  
  ##        posterior 75% quantile of k  
  ##        posterior mean of number of components M
  ##        posterior mode or medisn of number of components M
  ##        posterior 25% quantile of M  
  ##        posterior 75% quantile of M  
  ##        posterior probability of P(M_na | data) = 0
  
  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result.ARI = rep(0, MCMC_iteration-burn_in)
  result.cluster = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result.ARI[i-burn_in] <- adj.rand.index(History[[i]]$c, true_allocation)
    result.cluster[i-burn_in] <- History[[i]]$k
  }
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  
  Ad_matrix <- graph_from_adjacency_matrix(Ad_matrix, mode = "undirected")
  
  output <- c(modularity(Ad_matrix, result_c),
              mean(result.ARI),
              adj.rand.index(result_c, true_allocation),
              mean(result.cluster),
              median(result.cluster),
              as.numeric(quantile(result.cluster)[2]),
              as.numeric(quantile(result.cluster)[4]),
              mean(result_M),
              median(result_M),
              as.numeric(quantile(result_M)[2]),
              as.numeric(quantile(result_M)[4]),
              length(result_M_na[result_M_na==0]) / ( MCMC_iteration-burn_in )
  )
  
  return(output)
}

MCMC_output_dol <- function(MCMC_iteration, burn_in, History){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  
  ################################################################
  
  ## Return: A vector consisting of the following elements:
  ##        posterior mean of number of clusters k
  ##        posterior mode or medisn of number of clusters k
  ##        posterior 25% quantile of k  
  ##        posterior 75% quantile of k  
  ##        posterior mean of number of components M
  ##        posterior mode or medisn of number of components M
  ##        posterior 25% quantile of M  
  ##        posterior 75% quantile of M  
  ##        posterior probability of P(M_na | data) = 0
  
  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result.ARI = rep(0, MCMC_iteration-burn_in)
  result.cluster = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result.cluster[i-burn_in] <- History[[i]]$k
  }
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  
  output <- c(mean(result.cluster),
              median(result.cluster),
              as.numeric(quantile(result.cluster, probs = c(0.05))),
              as.numeric(quantile(result.cluster, probs = c(0.95))),
              mean(result_M),
              median(result_M),
              as.numeric(quantile(result_M, probs = c(0.05))),
              as.numeric(quantile(result_M, probs = c(0.95))),
              length(result_M_na[result_M_na==0]) / ( MCMC_iteration-burn_in )
  )
  
  return(output)
}

PP_output <- function(MCMC_iteration, burn_in, History, len, S){
  ## Input:
  ##        MCMC_iteration : Number of MCMC iterations
  ##        burn_in : burn-in period (number of iterations to discard).
  ##        History : A list used to store the results of the MCMC sampler, ex: History <- vector("list", MCMC_iteration) 
  ##        len : the truncation threshold: P(S \leq len : data)
  ##        S : S is a character string specifying either "M" or "k".
  ##          When S = "M", the function outputs the posterior probabilities of number of components M.
  ##          When S = "k", the function outputs the posterior probabilities of number of clusters k.
  ################################################################
  
  ## Return: A matrix storing of posterior probability of S
  
  result_prob <- rep(0, len)
  result <- rep(0, MCMC_iteration-burn_in)
  
  if(S=="M"){
    for (i in (burn_in+1):MCMC_iteration) {
      result[i-burn_in] <- History[[i]]$M
    } 
  }else if(S=="k"){
    for (i in (burn_in+1):MCMC_iteration) {
      result[i-burn_in] <- History[[i]]$k
    }
  }
  
  for (i in 1:len) {
    if(i < len){
      result_prob[i] <- length(result[result == i]) / (MCMC_iteration-burn_in)
    }else{
      result_prob[i] <- length(result[result >= len]) / (MCMC_iteration-burn_in)
    }
  }
  return(result_prob)
}

density_estimate <- function(grid_y, History, MCMC_iteration, burn_in){
  density_matrix <- matrix(0.0, length(grid_y), MCMC_iteration - burn_in)
  for (i in 1:length(grid_y)) {
    for (ii in (burn_in+1):MCMC_iteration) {
      if(History[[ii]]$M_na > 0){
        S_vector <- append(History[[ii]]$S_a, History[[ii]]$S_na)
        Theta_vector <- append(History[[ii]]$Theta_a, History[[ii]]$Theta_na)
        Sigma_vector <- append(History[[ii]]$Sigma_a, History[[ii]]$Sigma_na)
      }else{
        S_vector <- History[[ii]]$S_a
        Theta_vector <- History[[ii]]$Theta_a
        Sigma_vector <- History[[ii]]$Sigma_a
      }
      S_vector <- S_vector / sum(S_vector)
      density_matrix[i, ii-(burn_in)] <- sum(S_vector * dnorm(grid_y[i], Theta_vector, sqrt(Sigma_vector)))
    }
  }
  
  density_estimate <- data.frame(
    x = grid_y,
    posterior_mean = rep(0, length(grid_y)),
    quontile_0025 = rep(0, length(grid_y)),
    quontile_0975 = rep(0, length(grid_y))
  )
  
  density_estimate$posterior_mean <- apply(density_matrix, 1, mean)
  for (i in 1:length(grid_y)) {
    density_estimate$quontile_0025[i] <- quantile(density_matrix[i,], 0.025)
    density_estimate$quontile_0975[i] <- quantile(density_matrix[i,], 0.975)
  }
  
  return(density_estimate)
}
