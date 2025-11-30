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

sample_M_telescope_ig <- function(M_max, k, u, alp, lam, c){
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

MCMC_IG_network <- function(Ad_matrix, MCMC_iteration, Q_shape1, Q_shape2, a_lambda, b_lambda, Alpha){
  History <- vector("list", MCMC_iteration)
  node_size <- length(Ad_matrix[1,])
  data1 <- Ad_matrix
  data1[lower.tri(data1)] <- 0
  
  ### initial values ###
  # kk <- cluster_louvain(graph_from_adjacency_matrix(data, mode = "undirected"))
  # k <- length(unique(kk$membership))
  # M_na <- 0
  # M <- k
  # c <- kk$membership
  
  k <- 10
  M_na <- 0
  M <- k
  c <- rep(c(1:k), each = (node_size %/% k)+1)[1:node_size]
  
  S_a <- rep(1/M, k)
  Q <- matrix(0, M, M)
  for (r in 1:M) {
    for (s in r:M) {
      Q[r, s] <- rbeta(1, Q_shape1, Q_shape2)
      Q[s, r] <- Q[r, s]
    }
  }
  lambda = 1.0
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a,
                       Q = Q,
                       lambda = lambda)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector <- append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
    }else{
      S_vector <- History[[iter-1]]$S_a
    }
    
    ### update U_n ###
    U_n <- rgamma(1, shape=node_size, rate=sum(S_vector))
    
    ### update c ###
    c_old <- History[[iter-1]]$c
    M_old <- History[[iter-1]]$M
    Q_old <- History[[iter-1]]$Q
    c_update <- c_old
    for (i in 1:node_size) {
      weight_i <- rep(0, M_old)
      dummy_c <- c_update
      for (kk in 1:M_old) {
        dummy_c[i] <- kk
        weight_i[kk] <- S_vector[kk]*exp(loglike(dummy_c, Q_old, Ad_matrix, i, node_size))
        # exp(loglike(clusterAssign_temp,Q,data,i,n)
      }
      if(sum(weight_i) == 0){
        for (kk in 1:M_old) {
          weight_i[kk] <- log(S_vector[kk]) + loglike(Ad_matrix, dummy_c, Q_old, M_old, i)
        }
        weight_i <- exp(weight_i - max(weight_i))
      }
      c_update[i] <- sample(1:M_old, 1, replace = TRUE, prob = weight_i)
    }
    k <- length(unique(c_update))
    
    ### rename ###
    c_rename <- c_update
    c_unique <- sort(unique(c_update))
    S_a <- S_vector[c_unique]
    Q <- as.matrix(Q_old[c_unique, c_unique])
    for (i in 1:k) {
      c_rename[c_update == c_unique[i]] <- i
    }
    
    ### update lambda ###
    psi <- exp(log_psi(U_n, Alpha))
    a_ast <- k + a_lambda
    weights_lambda <- c(psi*(k + a_lambda - 1.0), k*(b_lambda + 1.0 - psi))
    randum_lambda <- sample(1:2, 1, replace = TRUE, prob = weights_lambda)
    if(randum_lambda == 1){
      lambda <- rgamma(1, shape = a_ast, rate = 1.0 - psi + b_lambda)
    }else{
      lambda <- rgamma(1, shape = a_ast - 1.0, rate = 1.0 - psi + b_lambda)
    }
    
    ### update M_na ###
    weight_M_na <- c(exp(log_psi(U_n, Alpha))*lambda, k)
    randum_M_na <- sample(1:2, 1, replace = TRUE, prob = weight_M_na)
    if(randum_M_na == 1){
      M_na <- rpois(1, exp(log_psi(U_n, Alpha))*lambda) + 1
    }else{
      M_na <- rpois(1, exp(log_psi(U_n, Alpha))*lambda)
    }
    
    ### update allocated parameters ###
    c_table <- table(c_rename)
    for (i in 1:k) {
      S_a[i] <- rgig(1, lambda = c_table[i] - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
    }
    
    AA <- matrix(0, k, k)
    NN <- matrix(0, k, k)
    for (r in 1:k){
      for (s in r:k)
      {
        AA[r,s] <- sum(data1[c_rename==r,c_rename==s]) + sum(data1[c_rename==s, c_rename==r]) - (r==s)*sum(data1[c_rename==s, c_rename==r])
        med <- matrix(0, node_size, node_size)
        med[which(c_rename==r),which(c_rename==s)] <- 1
        med1 <- matrix(0,node_size, node_size)
        med1[which(c_rename==s), which(c_rename==r)] <- 1
        NN[r,s] <- sum(med*lower.tri(med)) + sum(med1*lower.tri(med1))-(r==s)*sum(med1*lower.tri(med1))
        Q[r,s] <- rbeta(1, AA[r,s]+Q_shape1, NN[r,s]-AA[r,s]+Q_shape2)
        Q[s,r] <- Q[r,s]
      }
    }
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na <- rgig(M_na, lambda = - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
      QQ <- matrix(NA, k + M_na, k + M_na)
      QQ[1:k, 1:k] <- Q
      for (r in  (k+1):(k+M_na)) {
        for (s in 1:(k+M_na)) {
          QQ[r, s] <- rbeta(1, Q_shape1, Q_shape2)
          QQ[s, r] <- QQ[r, s]
        }
      }
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a, S_na=S_na,
                             Q = QQ,
                             lambda=lambda)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a,
                             Q=Q,
                             lambda=lambda)
    }
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

MCMC_IG_network_tele <- function(Ad_matrix, MCMC_iteration, Q_shape1, Q_shape2, a_lambda, b_lambda, Alpha, M_max){
  History <- vector("list", MCMC_iteration)
  node_size <- length(Ad_matrix[1,])
  data1 <- Ad_matrix
  data1[lower.tri(data1)] <- 0
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
  # kk <- cluster_louvain(graph_from_adjacency_matrix(data, mode = "undirected"))
  # k <- length(unique(kk$membership))
  # M_na <- 0
  # M <- k
  # c <- kk$membership
  
  k <- 10
  M_na <- 0
  M <- k
  c <- rep(c(1:k), each = (node_size %/% k)+1)[1:node_size]
  S_a <- rep(1/M, k)
  Q <- matrix(0, M, M)
  for (r in 1:M) {
    for (s in r:M) {
      Q[r, s] <- rbeta(1, Q_shape1, Q_shape2)
      Q[s, r] <- Q[r, s]
    }
  }
  lambda = 1.0
  Alpha_ini <- Alpha
  Alp <- Alpha
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a,
                       Q = Q,
                       lambda = lambda)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector <- append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
    }else{
      S_vector <- History[[iter-1]]$S_a
    }
    
    ### update U_n ###
    U_n <- rgamma(1, shape=node_size, rate=sum(S_vector))
    
    ### update c ###
    c_old <- History[[iter-1]]$c
    M_old <- History[[iter-1]]$M
    Q_old <- History[[iter-1]]$Q
    c_update <- c_old
    for (i in 1:node_size) {
      weight_i <- rep(0, M_old)
      dummy_c <- c_update
      for (kk in 1:M_old) {
        dummy_c[i] <- kk
        weight_i[kk] <- S_vector[kk]*exp(loglike(dummy_c, Q_old, Ad_matrix, i, node_size))
        # exp(loglike(clusterAssign_temp,Q,data,i,n)
      }
      if(sum(weight_i) == 0){
        for (kk in 1:M_old) {
          weight_i[kk] <- log(S_vector[kk]) + loglike(Ad_matrix, dummy_c, Q_old, M_old, i)
        }
        weight_i <- exp(weight_i - max(weight_i))
      }
      c_update[i] <- sample(1:M_old, 1, replace = TRUE, prob = weight_i)
    }
    k <- length(unique(c_update))
    
    ### rename ###
    c_rename <- c_update
    c_unique <- sort(unique(c_update))
    S_a <- S_vector[c_unique]
    Q <- as.matrix(Q_old[c_unique, c_unique])
    for (i in 1:k) {
      c_rename[c_update == c_unique[i]] <- i
    }
    
    ### update lambda ###
    psi <- exp(log_psi(U_n, Alpha))
    a_ast <- k + a_lambda
    weights_lambda <- c(psi*(k + a_lambda - 1.0), k*(b_lambda + 1.0 - psi))
    randum_lambda <- sample(1:2, 1, replace = TRUE, prob = weights_lambda)
    if(randum_lambda == 1){
      lambda <- rgamma(1, shape = a_ast, rate = 1.0 - psi + b_lambda)
    }else{
      lambda <- rgamma(1, shape = a_ast - 1.0, rate = 1.0 - psi + b_lambda)
    }
    
    ### update M and M_na ###
    M <- sample_M_telescope_ig(M_max, k, U_n, Alp, lambda, c_rename)
    M_na <- M - k
    
    MH_result_Alp <- sample_Alp_telescope(Alp, df1, df2, U_n, k, M, c_rename, s_pro)
    Alp <- MH_result_Alp$Alp
    acc <- MH_result_Alp$acc
    Alpha <- Alp / M
    
    ### update allocated parameters ###
    c_table <- table(c_rename)
    for (i in 1:k) {
      S_a[i] <- rgig(1, lambda = c_table[i] - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
    }
    
    AA <- matrix(0, k, k)
    NN <- matrix(0, k, k)
    for (r in 1:k){
      for (s in r:k)
      {
        AA[r,s] <- sum(data1[c_rename==r,c_rename==s]) + sum(data1[c_rename==s, c_rename==r]) - (r==s)*sum(data1[c_rename==s, c_rename==r])
        med <- matrix(0, node_size, node_size)
        med[which(c_rename==r),which(c_rename==s)] <- 1
        med1 <- matrix(0,node_size, node_size)
        med1[which(c_rename==s), which(c_rename==r)] <- 1
        NN[r,s] <- sum(med*lower.tri(med)) + sum(med1*lower.tri(med1))-(r==s)*sum(med1*lower.tri(med1))
        Q[r,s] <- rbeta(1, AA[r,s]+Q_shape1, NN[r,s]-AA[r,s]+Q_shape2)
        Q[s,r] <- Q[r,s]
      }
    }
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na <- rgig(M_na, lambda = - 0.5, chi = Alpha^2, psi = 2 * U_n + 1)
      QQ <- matrix(NA, k + M_na, k + M_na)
      QQ[1:k, 1:k] <- Q
      for (r in  (k+1):(k+M_na)) {
        for (s in 1:(k+M_na)) {
          QQ[r, s] <- rbeta(1, Q_shape1, Q_shape2)
          QQ[s, r] <- QQ[r, s]
        }
      }
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a, S_na=S_na,
                              Q = QQ,
                              lambda=lambda,
                              shape = Alpha,
                              acc = acc)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Q=Q,
                              lambda=lambda,
                              shape = Alpha,
                              acc = acc)
    }
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

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

MCMC_output <- function(MCMC_iteration, burn_in, History, true_allocation, true_M){
  result_M <- rep(0, MCMC_iteration-burn_in)
  result_M_na <- rep(0, MCMC_iteration-burn_in)
  result_k <- rep(0, MCMC_iteration-burn_in)
  result_lambda <- rep(0, MCMC_iteration-burn_in)
  result_post_ari <- rep(0, MCMC_iteration-burn_in)
  result_shape <- rep(0, MCMC_iteration-burn_in)
  
  
  for (i in (burn_in+1):MCMC_iteration) {
    result_M[i-burn_in] <- History[[i]]$M
    result_M_na[i-burn_in] <- History[[i]]$M_na
    result_k[i-burn_in] <- History[[i]]$k
    result_lambda[i-burn_in] <- History[[i]]$lambda
    result_post_ari[i-burn_in] <- adj.rand.index(History[[i]]$c, true_allocation)
    result_shape[i-burn_in] <- History[[i]]$shape
  }
  
  fit.MFM <- list(Iterates = History)
  result.MFM <- getDahl(fit.MFM, burn = burn_in)
  result_c <- result.MFM$c
  output <- c(mean(result_M),
              length(result_M_na[result_M_na==0]) / (MCMC_iteration-burn_in),
              mean(result_k),
              mean(result_post_ari),
              adj.rand.index(result_c, true_allocation),
              mean(result_lambda),
              mean(result_shape))
  
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