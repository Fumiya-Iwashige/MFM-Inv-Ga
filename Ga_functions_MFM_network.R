BG_Ga_network_MFM_Poi <- function(Ad_matrix, MCMC_iteration, Q_shape1, Q_shape2, lambda){
  ## Model: M-1 \sim Poisson(lambda)
  ##        A_ij | c,Q \sim Ber(Q_{c_i,c_j})
  ##        Q_{rs} \sim Beta(Q_shape1, Q_shape2), r,s = 1,...,k
  ##        h_m | gamma \sim Gamma(gamma, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
  ##        
  
  ################################################################
  
  ## Input: Ad_matrix : An adjacent matrix of size (node size × node size).
  ##        MCMC_iteration : Number of MCMC iterations
  ##        Q_shape1, Q_shape2 : hyperparameters for Beta prior of Q
  ##        lambda : hyperparameter for the Poisson prior of M-1
  ##        k_ini : initial number of occupied clusters.
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        k : number of clusters
  ##        c : label vector
  ##        Q : probability matrix
  ##        S_a, S_na: un-normalized component weights
  
  History <- vector("list", MCMC_iteration)
  node_size <- length(Ad_matrix[1,])
  data1 <- Ad_matrix
  data1[lower.tri(data1)] <- 0
  
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
  gamma <- 1.0
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
    
    ### update M_na ###
    weight_M_na <- c(exp(log_psi_ga(U_n, gamma))*lambda, k)
    randum_M_na <- sample(1:2, 1, replace = TRUE, prob = weight_M_na)
    if(randum_M_na == 1){
      M_na <- rpois(1, exp(log_psi_ga(U_n, gamma))*lambda) + 1
    }else{
      M_na <- rpois(1, exp(log_psi_ga(U_n, gamma))*lambda)
    }
    M <- k + M_na
    
    ### update shape gamma ###
    MH_result_Ga <- sample_Ga_sta(gamma, df1, df2, U_n, k, M, c_rename, s_pro)
    gamma <- MH_result_Ga$Ga
    acc <- MH_result_Ga$acc
    
    ### update allocated parameters ###
    c_table <- table(c_rename)
    for (i in 1:k) {
      S_a[i] <- rgamma(1, shape = c_table[i]+gamma, rate=U_n+1)
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
      S_na <- rgamma(M_na, shape = gamma, rate = U_n+1)
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
                              Q = QQ)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Q=Q)
    }
    # cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

BG_Ga_network_MFM_NB <- function(Ad_matrix, MCMC_iteration, Q_shape1, Q_shape2, a_lambda, b_lambda){
  ## Model: M-1 | lambda \sim Poisson(lambda)
  ##        lambda \sim Gamma(a_lambda, b_lambda)
  ##        A_ij | c,Q \sim Ber(Q_{c_i,c_j})
  ##        Q_{rs} \sim Beta(Q_shape1, Q_shape2), r,s = 1,...,k
  ##        h_m | gamma \sim Gamma(gamma, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
  ##        
  
  ################################################################
  
  ## Input: Ad_matrix : An adjacent matrix of size (node size × node size).
  ##        MCMC_iteration : Number of MCMC iterations
  ##        Q_shape1, Q_shape2 : hyperparameters for Beta prior of Q
  ##        a_lambda, b_lambda : hyperparameter for the Gamma prior of lambda
  ##        k_ini : initial number of occupied clusters.
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        k : number of clusters
  ##        c : label vector
  ##        Q : probability matrix
  ##        S_a, S_na: un-normalized component weights
  
  History <- vector("list", MCMC_iteration)
  node_size <- length(Ad_matrix[1,])
  data1 <- Ad_matrix
  data1[lower.tri(data1)] <- 0
  
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
  gamma <- 1.0
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
  lambda <- 1.0
  
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
    psi <- exp(log_psi_ga(U_n, gamma))
    a_ast <- k + a_lambda
    weights_lambda <- c(psi*(k + a_lambda - 1.0), k*(b_lambda + 1.0 - psi))
    randum_lambda <- sample(1:2, 1, replace = TRUE, prob = weights_lambda)
    if(randum_lambda == 1){
      lambda <- rgamma(1, shape = a_ast, rate = 1.0 - psi + b_lambda)
    }else{
      lambda <- rgamma(1, shape = a_ast - 1.0, rate = 1.0 - psi + b_lambda)
    }
    
    ### update M_na ###
    weight_M_na <- c(exp(log_psi_ga(U_n, gamma))*lambda, k)
    randum_M_na <- sample(1:2, 1, replace = TRUE, prob = weight_M_na)
    if(randum_M_na == 1){
      M_na <- rpois(1, exp(log_psi_ga(U_n, gamma))*lambda) + 1
    }else{
      M_na <- rpois(1, exp(log_psi_ga(U_n, gamma))*lambda)
    }
    M <- k + M_na
    
    ### update shape gamma ###
    MH_result_Ga <- sample_Ga_sta(gamma, df1, df2, U_n, k, M, c_rename, s_pro)
    gamma <- MH_result_Ga$Ga
    acc <- MH_result_Ga$acc
    
    ### update allocated parameters ###
    c_table <- table(c_rename)
    for (i in 1:k) {
      S_a[i] <- rgamma(1, shape = c_table[i]+gamma, rate=U_n+1)
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
      S_na <- rgamma(M_na, shape = gamma, rate = U_n+1)
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
                              Q = QQ)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Q=Q)
    }
    # cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}
