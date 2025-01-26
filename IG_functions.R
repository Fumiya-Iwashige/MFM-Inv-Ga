log_psi = function(u, Alpha) {Alpha * (1 - sqrt(1 + 2 * u))}

alpha_u = function(u){1 + 2 * u}

log_cummulants_ul = function(u, Alpha, l){
  result = 0
  x = Alpha*sqrt(alpha_u(u))
  if(x >= 100){
    result = l * log(Alpha) + log_psi(u, Alpha) - l * log(alpha_u(u, Alpha)) * 0.5 
  }else{
    result = l * log(Alpha) + log_psi(u, Alpha) - l * log(alpha_u(u)) * 0.5 + log(BesselK(l - 0.5, x)) - log(BesselK(0.5, x))
  }
  return(result)
}

log_Bessel_ratio = function(u, Alpha, l){
  result = 0
  x = Alpha*sqrt(alpha_u(u))
  if(BesselK(l - 0.5, x) > 0 & BesselK(0.5, x) > 0){
    result = log(BesselK(l - 0.5, x)) - log(BesselK(0.5, x))
  }else{
    result = 0
  }
  return(result)
}

Psi = function(u, Alpha, lambda, k){ lambda^(k -1) * (lambda * exp(log_psi(u, Alpha)) + k) * exp(lambda * (exp(log_psi(u, Alpha)) - 1)) }

MCMC_IG_blocking_multi_normal = function(data, MCMC_iteration, m_0, c_0, C_0, a_lambda, b_lambda, Alpha){
  History = vector("list", MCMC_iteration)
  data_size = length(data[1,])
  dimension = length(data[,1])

  ### initial values ###
  K = 5
  k = sample(1:K, 1, replace = TRUE, prob = rep(1/K, K))
  M_na = sample(1:K, 1, replace = TRUE, prob = rep(1/K, K))
  M = k + M_na
  c = sample(1:k, data_size, replace = TRUE, prob = rep(1/k, k)); c[1:k] = 1:k
  S_a = runif(k, 0, 1)
  Theta_a = matrix(0, dimension, k)
  Sigma_a = array(0, dim = c(dimension, dimension, k))
  for (i in 1:k) {
    Sigma_a[,,i] = rinvwishart(c_0, C_0)
    Theta_a[,i] = rmvnorm(1, m_0, Sigma_a[,,i])
  }
  S_na = runif(M_na, 0, 1)
  Theta_na = matrix(0, dimension, M_na)
  Sigma_na = array(0, dim = c(dimension, dimension, M_na))
  for (i in 1:M_na) {
    Sigma_na[,,i] = rinvwishart(c_0, C_0)
    Theta_na[,i] = rmvnorm(1, m_0, Sigma_na[,,i])
  }
  lambda = rgamma(1, a_lambda, b_lambda)
  X = matrix(1, 1, data_size)
  
  if(Alpha < 1e-2){
    S_a = rinvgauss(k, Alpha, Alpha^2)
    S_na = rinvgauss(M_na, Alpha, Alpha^2)
  }
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, S_na=S_na,
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       Theta_na=Theta_na, Sigma_na=Sigma_na,
                       lambda = lambda)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector = append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
      Theta_matrix = abind(History[[iter-1]]$Theta_a, History[[iter-1]]$Theta_na)
      Sigma_array = abind(History[[iter-1]]$Sigma_a, History[[iter-1]]$Sigma_na)
    }else{
      S_vector = History[[iter-1]]$S_a
      Theta_matrix = History[[iter-1]]$Theta_a
      Sigma_array = History[[iter-1]]$Sigma_a
    }
    
    ### update U_n ###
    U_n = rgamma(1, shape=data_size, rate=sum(S_vector))
    
    ### update c ###
    c_old = History[[iter-1]]$c
    M_old = History[[iter-1]]$M
    c_update = rep(0, data_size)
    sample_vector = sample(1:data_size)
    for (i in sample_vector) {
      weight_i = rep(0, M_old)
      for (j in 1:M_old) {
        weight_i[j] = S_vector[j]*dmvnorm(data[,i], Theta_matrix[,j], Sigma_array[,,j], log=FALSE)
      }
      if(sum(weight_i) == 0){
        for (j in 1:M_old) {
          weight_i[j] = log(S_vector[j]) + dmvnorm(data[,i], Theta_matrix[,j], Sigma_array[,,j], log=TRUE)
        }
        weight_i = exp(weight_i - max(weight_i))
      }
      c_update[i] = sample(1:M_old, 1, replace = TRUE, prob = weight_i)
    }
    k = length(unique(c_update))
    
    ### rename ###
    c_rename = c_update
    c_unique = sort(unique(c_update))
    S_a = S_vector[c_unique]
    Theta_a = as.matrix(Theta_matrix[,c_unique])
    Sigma_a = array(Sigma_array[,,c_unique], dim = c(dimension, dimension, k))
    for (i in 1:k) {
      c_rename[c_update == c_unique[i]] = i
    }
    
    ### update lambda ###
    psi = exp(log_psi(U_n, Alpha))
    a_ast = k + a_lambda
    weights_lambda = c(psi*(k + a_lambda - 1.0), k*(b_lambda + 1.0 - psi))
    randum_lambda = sample(1:2, 1, replace = TRUE, prob = weights_lambda)
    if(randum_lambda == 1){
      lambda = rgamma(1, shape = a_ast, rate = 1.0 - psi + b_lambda)
    }else{
      lambda = rgamma(1, shape = a_ast - 1.0, rate = 1.0 - psi + b_lambda)
    }
    
    ### update M_na ###
    weight_M_na = c(exp(log_psi(U_n, Alpha))*lambda, k)
    randum_M_na = sample(1:2, 1, replace = TRUE, prob = weight_M_na)
    if(randum_M_na == 1){
      M_na = rpois(1, exp(log_psi(U_n, Alpha))*lambda) + 1
    }else{
      M_na = rpois(1, exp(log_psi(U_n, Alpha))*lambda)
    }
    
    ### update allocated parameters:S_a, Theta_a and Sigma_a ###
    c_table = table(c_rename)
    sample_vector = sample(1:k)
    for (i in sample_vector) {
      S_a[i] = rgig(1, lambda = c_table[i] - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
      
      N_i = c_table[i]
      Y_i = data[, which(c_rename == i)]
      X_i = X[which(c_rename == i)]
      if(length(X_i) == 1){
        X_i = as.matrix( X[which(c_rename == i)])
        Y_i = as.matrix(Y_i, 1, 2)
      }
      B_i = 1 / (N_i + 1)
      M_i = (Y_i %*% X_i + m_0) * B_i
      c_i = c_0 + N_i
      
      Theta_a[,i] = rmvnorm(1, M_i, Sigma_a[,,i]*B_i)
      Sigma_a[,,i] = rinvwishart(c_i, C_0 + tcrossprod(Y_i - M_i %*% t(X_i)) + tcrossprod(M_i - m_0))
    }
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na = rgig(M_na, lambda = - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
      Theta_na = matrix(0,dimension, M_na)
      Sigma_na = array(0, dim = c(dimension, dimension, M_na))
      sample_vector = sample(1:M_na)
      for (i in sample_vector) {
        Sigma_na[,,i] = rinvwishart(c_0, C_0)
        Theta_na[,i] = rmvnorm(1, m_0, Sigma_na[,,i])
      }
      
      History[[iter]] = list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a, S_na=S_na,
                             Theta_a=Theta_a, Theta_na=Theta_na,
                             Sigma_a=Sigma_a, Sigma_na=Sigma_na,
                             lambda=lambda)
    }else{
      History[[iter]] = list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a,
                             Theta_a=Theta_a,
                             Sigma_a=Sigma_a,
                             lambda=lambda)
    }
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

MCMC_IG_blocking_uni_normal = function(data, MCMC_iteration, m_0, c_0, d_0, D_0, w_0, W_0, a_lambda, b_lambda, Alpha){
  History = vector("list", MCMC_iteration)
  data_size = length(data)
  m = Alpha
  xi = Alpha^2
  
  ### initial values ###
  K = 20
  k = sample(1:K, 1, replace = TRUE, prob = rep(1/K, K))
  M_na = sample(1:K, 1, replace = TRUE, prob = rep(1/K, K))
  M = k + M_na
  c = sample(1:k, data_size, replace = TRUE, prob = rep(1/k, k)); c[1:k] = 1:k 
  S_a = rinvgauss(k, Alpha, Alpha^2)
  Theta_a = rnorm(k, 0, 1)
  Sigma_a = rinvgamma(k, 1, 1)
  S_na = rinvgauss(M_na, Alpha, Alpha^2)
  Theta_na = rnorm(M_na, 0, 1)
  Sigma_na = rinvgamma(M_na, 1, 1)
  X = matrix(1, 1, data_size)
  lambda = rgamma(1, a_lambda, b_lambda)
  B_0 = rgamma(1, 1, 1)
  C_0 = rgamma(1, 1, 1)
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, S_na=S_na,
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       Theta_na=Theta_na, Sigma_na=Sigma_na,
                       lambda = lambda,
                       B_0 = B_0,
                       C_0 = C_0)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector = append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
      Theta_vector = append(History[[iter-1]]$Theta_a, History[[iter-1]]$Theta_na)
      Sigma_vector = append(History[[iter-1]]$Sigma_a, History[[iter-1]]$Sigma_na)
    }else{
      S_vector = History[[iter-1]]$S_a
      Theta_vector = History[[iter-1]]$Theta_a
      Sigma_vector = History[[iter-1]]$Sigma_a
    }
    
    ### update U_n ###
    U_n = rgamma(1, shape=data_size, rate=sum(S_vector))
    
    ### update c ###
    c_old = History[[iter-1]]$c
    M_old = History[[iter-1]]$M
    c_update = c_old
    sample_vector = sample(1:data_size)
    for (i in sample_vector) {
      weight_i = rep(0, M_old)
      for (j in 1:M_old) {
        weight_i[j] = S_vector[j] * dnorm(data[i], Theta_vector[j], sqrt(Sigma_vector[j]), log=FALSE)
      }
      if(sum(weight_i) == 0){
        for (j in 1:M_old) {
          weight_i[j] = log(S_vector[j]) + dnorm(data[i], Theta_vector[j], sqrt(Sigma_vector[j]), log=TRUE) + 2 * log(10)
        }
        weight_i = exp(weight_i - max(weight_i))
      }
      c_update[i] = sample(1:M_old, 1, replace = TRUE, prob = weight_i)
    }
    k = length(unique(c_update))
    
    ### rename ###
    c_rename = c_update
    c_unique = sort(unique(c_update))
    S_a = S_vector[c_unique]
    Theta_a = Theta_vector[c_unique]
    Sigma_a = Sigma_vector[c_unique]
    for (i in 1:k) {
      c_rename[c_update == c_unique[i]] = i
    }
    
    ### update lambda ###
    psi = exp(log_psi(U_n, Alpha))
    a_ast = k + a_lambda
    weights_lambda = c(psi*(k + a_lambda - 1.0), k*(b_lambda + 1.0 - psi))
    randum_lambda = sample(1:2, 1, replace = TRUE, prob = weights_lambda)
    if(randum_lambda == 1){
      lambda = rgamma(1, shape = a_ast, rate = 1.0 - psi + b_lambda)
    }else{
      lambda = rgamma(1, shape = a_ast - 1.0, rate = 1.0 - psi + b_lambda)
    }
    
    ### update M_na ###
    weight_M_na = c(exp(log_psi(U_n, Alpha))*lambda, k)
    randum_M_na = sample(1:2, 1, replace = TRUE, prob = weight_M_na)
    if(randum_M_na == 1){
      M_na = rpois(1, exp(log_psi(U_n, Alpha))*lambda) + 1
    }else{
      M_na = rpois(1, exp(log_psi(U_n, Alpha))*lambda)
    }
    
    ### update allocated parameters:S_a, Theta_a and Sigma_a ###
    c_table = table(c_rename)
    sample_vector = sample(1:k)
    for (i in sample_vector) {
      S_a[i] = rgig(1, lambda = c_table[i] - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
      
      N_i = c_table[i]
      Y_i = data[which(c_rename == i)]
      X_i = X[which(c_rename == i)]
      B_i = 1 / (N_i + B_0)
      M_i = (sum(Y_i) + B_0 * m_0) * B_i
      c_i = c_0 + N_i/2
      C_k = B_0 * (Theta_a[i] - m_0)^2 + sum((Y_i - M_i)^2)
      
      Theta_a[i] = rnorm(1, M_i, sqrt(Sigma_a[i]*B_i))
      Sigma_a[i] = rinvgamma(1, c_i, C_0 + C_k / 2)
    }
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na = rgig(M_na, lambda = -0.5, chi = Alpha^2, psi = 2 * U_n + 1)
      Theta_na = rep(0, M_na)
      Sigma_na = rep(0, M_na)
      sample_vector = sample(1:M_na)
      for (i in sample_vector) {
        Sigma_na[i] = rinvgamma(1, c_0, C_0)
        Theta_na[i] = rnorm(1, m_0, sqrt(Sigma_na[i]/B_0))
      }
      
      ### update B_0 and C_0 ###
      W_1 = (k+M_na) /2 + w_0
      W_2 = W_0 + sum((append(Theta_a, Theta_na) - m_0)^2 / append(Sigma_a, Sigma_na)) * 0.5
      B_0 = rgamma(1, W_1, W_2)
      
      D_1 = (k + M_na) * c_0 + d_0
      D_2 = sum(1/append(Sigma_a, Sigma_na)) + D_0
      C_0 = rgamma(1, D_1, D_2)
      
      History[[iter]] = list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a, S_na=S_na,
                             Theta_a=Theta_a, Theta_na=Theta_na,
                             Sigma_a=Sigma_a, Sigma_na=Sigma_na,
                             lambda=lambda,
                             B_0 = B_0,
                             C_0 = C_0)
    }else{
      ### update B_0 ###
      W_1 = k /2 + w_0
      W_2 = W_0 + sum((Theta_a - m_0)^2 / Sigma_a) * 0.5
      B_0 = rgamma(1, W_1, W_2)
      
      D_1 = k * c_0 + d_0
      D_2 = sum(1/Sigma_a) + D_0
      C_0 = rgamma(1, D_1, D_2)
      
      History[[iter]] = list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a,
                             Theta_a=Theta_a,
                             Sigma_a=Sigma_a,
                             lambda=lambda,
                             B_0 = B_0,
                             C_0 = C_0)
    }
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}


## Dahl's method to summarize the samples from the MCMC
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

MCMC_output = function(MCMC_iteration, burn_in, History, true_allocation, true_M){
  result_M = rep(0, MCMC_iteration-burn_in)
  result_M_na = rep(0, MCMC_iteration-burn_in)
  result_lambda = rep(0, MCMC_iteration-burn_in)
  result_U = rep(0, MCMC_iteration-burn_in)
  result_post_ri = rep(0, MCMC_iteration-burn_in)
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] = History[[i]]$M
    result_M_na[i-burn_in] = History[[i]]$M_na
    result_lambda[i-burn_in] = History[[i]]$lambda
    result_U[i-burn_in] = History[[i]]$U
    result_post_ri[i-burn_in] = rand.index(History[[i]]$c, true_allocation)
  }
  
  fit.MFM = list(Iterates = History)
  result.MFM = getDahl(fit.MFM, burn = burn_in)
  result_c = result.MFM$c
  
  output = c(mean(result_M),
             (mean(result_M) - true_M)^2,
             length(result_M_na[result_M_na==0]) / (MCMC_iteration-burn_in),
             mean(result_post_ri),
             rand.index(result_c, true_allocation),
             mean(result_lambda),
             mean(result_U))
  
  return(output)
}

PP_output = function(MCMC_iteration, burn_in, History, len, S){
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