log_psi <- function(u, Alpha) {Alpha * (1 - sqrt(1 + 2 * u))}

alpha_u <- function(u){1 + 2 * u}

log_cummulant_IG <- function(u, Alpha, n){
  result <- 0.0
  hu <- sqrt(1 + 2 * u)
  z <- Alpha * hu
  log_psi_IG <- log_psi(u, Alpha)
  
  ### the log ratio of K_n-1/2(z) / K_1/2(z) ###
  value_bessel <- Re(Bessel::BesselK(z = z, nu = n - 0.5, expon.scaled = TRUE))
  if(Re(value_bessel) == Inf){
    log_K_num <- Re(Bessel::besselK.nuAsym(x = z, nu = n-0.5, k.max = 5, expon.scaled = TRUE, log = TRUE)) - z
  }else{
    log_K_num <- log(value_bessel) - z
  }
  logK_den <- 0.5 * (log(pi) - log(2) - log(z)) - z
  
  log_ratio_Bessel <- log_K_num - logK_den
  
  result <- n * log(Alpha) + log_psi_IG - ( n/2 ) * log(1 + 2*u) + log_ratio_Bessel
  return(result)
}

Psi <- function(u, Alpha, lambda, k){ lambda^(k -1) * (lambda * exp(log_psi(u, Alpha)) + k) * exp(lambda * (exp(log_psi(u, Alpha)) - 1)) }

sample_Alp_sta <- function(Alp, df1, df2, u, k, m, c, s_pro){
  c_table <- tabulate(c)
  log_post_Alp <- function(x){
    log_cummulant_vec <- rep(0.0, k)
    for (i in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_IG(u, x, c_table[i])
    }
    df(x=x, df1, df2, log=TRUE) + sum(log_cummulant_vec) + (m-k)*log_psi(u, x)
  }
  log_Alp_pro <- log(Alp) + rnorm(1, 0, s_pro)
  Alp_pro <- exp(log_Alp_pro)
  
  log_accept <- log_post_Alp(Alp_pro) - log_post_Alp(Alp) + log(Alp_pro) - log(Alp)
  accept <- Re(exp(log_accept))
  acc <- FALSE
  if (runif(1) <= accept) {
    Alp <- Alp_pro
    acc <- TRUE
  }
  return(list(Alp = Alp, acc = acc))
}

sample_Alp_telescope <- function(Alp, df1, df2, u, k, m, c, s_pro){
  c_table <- tabulate(c)
  log_post_Alp <- function(x){
    log_cummulant_vec <- rep(0.0, k)
    for (i in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_IG(u, x/m, c_table[i])
    }
    df(x=x, df1, df2, log=TRUE) + sum(log_cummulant_vec) + (m-k)*log_psi(u, x/m)
  }
  log_Alp_pro <- log(Alp) + rnorm(1, 0, s_pro)
  Alp_pro <- exp(log_Alp_pro)
  
  log_accept <- log_post_Alp(Alp_pro) - log_post_Alp(Alp) + log(Alp_pro) - log(Alp)
  accept <- Re(exp(log_accept))
  acc <- FALSE
  if (runif(1) <= accept) {
    Alp <- Alp_pro
    acc <- TRUE
  }
  return(list(Alp = Alp, acc = acc))
}

log_post_blambda <- function(b_lambda, a_lambda, lambda, a_p, b_p){
  dgamma(lambda, shape = a_lambda, rate = b_lambda, log = TRUE) + dbetapr(b_lambda, shape1=a_p, shape2=b_p, scale = 1, log = TRUE)
}

sample_M_IG_BNB <- function(M_max, k, u, alp, al_lam, a_pi, b_pi, c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- tabulate(c)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi(u, alp) + dbnbinom(x-1, al_lam, alpha = a_pi, beta = b_pi, log = TRUE)
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max)
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_M_TS_IG_Poi <- function(M_max, k, u, alp, lam, c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- tabulate(c)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi(u, alp/x) + dpois(x-1, lam, log = TRUE)
  }
  
  log_cummulant_vec <- rep(0.0, M_max)
  for (i in k:M_max) {
    for (j in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_vec[i] + log_cummulant_IG(u, alp/i, c_table[j])
    }
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max) + log_cummulant_vec[k:M_max]
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_M_TS_IG_NB <- function(M_max, k, u, alp, r, p, c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- tabulate(c)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi(u, alp/x) + dnbinom(x-1, size=r, prob=1-p, log=TRUE)
  }
  
  log_cummulant_vec <- rep(0.0, M_max)
  for (i in k:M_max) {
    for (j in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_vec[i] + log_cummulant_IG(u, alp/i, c_table[j])
    }
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max) + log_cummulant_vec[k:M_max]
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_M_TS_IG_BNB <- function(M_max, k, u, alp, al_lam, a_pi, b_pi,c){
  logprob_vec <- rep(-Inf, M_max)
  c_table <- tabulate(c)
  lquot <- function(x) {
    lgamma(x + 1) - lgamma(x - k + 1) + (x-k) * log_psi(u, alp/x) + dbnbinom(x-1, al_lam, alpha = a_pi, beta = b_pi, log = TRUE)
  }
  
  log_cummulant_vec <- rep(0.0, M_max)
  for (i in k:M_max) {
    for (j in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_vec[i] + log_cummulant_IG(u, alp/i, c_table[j])
    }
  }
  
  logprob_vec[k:M_max] <- lquot(k:M_max) + log_cummulant_vec[k:M_max]
  logprob_vec[k:M_max]  <- logprob_vec[k:M_max] - max(logprob_vec[k:M_max])
  prob_vec <- exp(logprob_vec)
  result <- sample(1:M_max, 1, prob = prob_vec, replace = TRUE)
  return(result)
}

sample_Alp_telescope <- function(Alp, df1, df2, u, k, m, c, s_pro){
  c_table <- tabulate(c)
  log_post_Alp <- function(x){
    log_cummulant_vec <- rep(0.0, k)
    for (i in 1:k) {
      log_cummulant_vec[i] <- log_cummulant_IG(u, x/m, c_table[i])
    }
    df(x=x, df1, df2, log=TRUE) + sum(log_cummulant_vec) + (m-k)*log_psi(u, x/m)
  }
  log_Alp_pro <- log(Alp) + rnorm(1, 0, s_pro)
  Alp_pro <- exp(log_Alp_pro)
  
  log_accept <- log_post_Alp(Alp_pro) - log_post_Alp(Alp) + log(Alp_pro) - log(Alp)
  accept <- Re(exp(log_accept))
  acc <- FALSE
  if (runif(1) <= accept) {
    Alp <- Alp_pro
    acc <- TRUE
  }
  return(list(Alp = Alp, acc = acc))
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