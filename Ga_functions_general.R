log_psi_ga <- function(u, gamma){ -gamma*log(1.0 + u) }
log_cummulants_Ga <- function(u, gamma, n){ -(n + gamma)*log(u + 1) + lgamma(gamma + n) - lgamma(gamma) }
Psi <- function(u, gamma, lambda, k){ lambda^(k -1) * (lambda * exp(log_psi_ga(u, gamma)) + k) * exp(lambda * (exp(log_psi_ga(u, gamma)) - 1)) }

sample_Ga_sta <- function(Ga, df1, df2, u, k, m, c, s_pro){
  c_table <- table(c)
  log_post_Ga <- function(x){
    df(x=x, df1, df2, log=TRUE) + sum(log_cummulants_Ga(u, x, c_table)) + (m-k)*log_psi_ga(u, x)
  }
  log_Ga_pro <- log(Ga) + rnorm(1, 0, s_pro)
  Ga_pro <- exp(log_Ga_pro)
  
  log_accept <- log_post_Ga(Ga_pro) - log_post_Ga(Ga) + log(Ga_pro) - log(Ga)
  accept <- exp(log_accept)
  acc <- FALSE
  if (runif(1) <= accept) {
    Ga <- Ga_pro
    acc <- TRUE
  }
  return(list(Ga = Ga, acc = acc))
}

sample_Ga_telescope <- function(Ga, df1, df2, u, k, m, c, s_pro){
  c_table <- table(c)
  log_post_Ga <- function(x){
    df(x=x, df1, df2, log=TRUE) + sum(log_cummulants_Ga(u, x/m, c_table)) + (m-k)*log_psi_ga(u, x/m)
  }
  log_Ga_pro <- log(Ga) + rnorm(1, 0, s_pro)
  Ga_pro <- exp(log_Ga_pro)
  
  log_accept <- log_post_Ga(Ga_pro) - log_post_Ga(Ga) + log(Ga_pro) - log(Ga)
  accept <- exp(log_accept)
  acc <- FALSE
  if (runif(1) <= accept) {
    Ga <- Ga_pro
    acc <- TRUE
  }
  return(list(Ga = Ga, acc = acc))
}

sample_M_Ga_BNB <- function(M_max, k, u, ga, al_lam, a_pi, b_pi,c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- table(c)
  log_psi_u <- log_psi_ga(u, ga)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi_ga(u, ga) + dbnbinom(x-1, al_lam, alpha = a_pi, beta = b_pi, log = TRUE)
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max)
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_M_TS_Ga_Poi <- function(M_max, k, u, ga, lam, c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- table(c)
  log_psi_u <- log_psi_ga(u, ga)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi_ga(u, ga/x) + dpois(x-1, lam, log = TRUE)
  }
  
  log_cummulant_vec <- rep(0.0, M_max)
  for (i in k:M_max) {
    for (j in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_vec[i] + log_cummulants_Ga(u, ga/i, c_table[j])
    }
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max) + log_cummulant_vec[k:M_max]
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_M_TS_Ga_NB <- function(M_max, k, u, ga, r, p, c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- table(c)
  #log_psi_u <- log_psi_ga(u, ga)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi_ga(u, ga/x) + dnbinom(x-1, size=r, prob=1-p, log=TRUE)
  }
  
  log_cummulant_vec <- rep(0.0, M_max)
  for (i in k:M_max) {
    for (j in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_vec[i] + log_cummulants_Ga(u, ga/i, c_table[j])
    }
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max) + log_cummulant_vec[k:M_max]
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_M_TS_Ga_BNB <- function(M_max, k, u, ga, al_lam, a_pi, b_pi,c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- table(c)
  log_psi_u <- log_psi_ga(u, ga)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi_ga(u, ga/x) + dbnbinom(x-1, al_lam, alpha = a_pi, beta = b_pi, log = TRUE)
  }
  
  log_cummulant_vec <- rep(0.0, M_max)
  for (i in k:M_max) {
    for (j in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_vec[i] + log_cummulants_Ga(u, ga/i, c_table[j])
    }
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max) + log_cummulant_vec[k:M_max]
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

log_post_blambda <- function(b_lambda, a_lambda, lambda, a_p, b_p){
  dgamma(lambda, shape = a_lambda, rate = b_lambda, log = TRUE) + dbetapr(b_lambda, shape1=a_p, shape2=b_p, scale = 1, log = TRUE)
}

loglike <- function(clusterassign,param,data,J,n)
{
  clustersize = max(clusterassign)
  param = as.matrix(param)
  
  if (J==1) {result2 = 0
  for (ii in c((J+1):n))
  {
    result2 = result2 + data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
  }
  output = sum(result2)} else if (J==n){
    result = 0
    for (ii in c(1:(J-1)))
    {
      result = result + data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    output = sum(result)
  } else {
    result = 0
    for (ii in c(1:(J-1)))
    {
      result = result + data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    
    result2 = 0
    for (ii in c((J+1):n))
      
    {
      result2 = result2 + data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
    }
    output = sum(result)+sum(result2)}
  output
}