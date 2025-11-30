TS_Ga_MultiNormal_DMFM_Poi <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, lambda, k_ini, M_max){
  ## Model: M-1 \sim Poisson(lambda)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma, M \sim Gamma(gamma/M, 1), m=1,...,M
  ##        Alpha \sim F(df1, df2)
  ##        
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        b_0 : prior mean vector of N_r(b_0, B_0)
  ##        B_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        c_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        g_0, G_0 : hyperparameter for the Wishart prior of C.
  ##        lambda : hyperparameter for the Poisson prior of M-1
  ##        k_ini : initial number of occupied clusters.
  ##        M_max : upper bound for sampling M
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means, r × k and r × M_na
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        C_0 : hierarchical positive definite symmetric matrix
  ##        shape : shape parameter of h
  History <- vector("list", MCMC_iteration)
  data_size <- length(data[1,])
  dimension <- length(data[,1])
  B_0_inv <- chol2inv(chol(B_0))
  
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
  k <- k_ini
  M_na <- 0
  M <- k + M_na
  gamma <- 1.0
  Ga <- 1.0
  c_kmeans <- kmeans(t(data), centers = k, nstart = 100)
  c <- c_kmeans$cluster
  
  C_0 <- g_0 * chol2inv(chol(G_0))
  C_0_inv <- chol2inv(chol(C_0))
  
  S_a <- rep(1/M, k)
  Theta_a <- t(c_kmeans$centers)
  Sigma_a <- array(0, dim = c(dimension, dimension, k))
  Sigma_a[, , 1:k] <- 0.5 * C_0
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, 
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       C_0 = C_0,
                       lambda = lambda,
                       shape = gamma)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector <- append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
      Theta_matrix <- abind(History[[iter-1]]$Theta_a, History[[iter-1]]$Theta_na)
      Sigma_array <- abind(History[[iter-1]]$Sigma_a, History[[iter-1]]$Sigma_na)
    }else{
      S_vector <- History[[iter-1]]$S_a
      Theta_matrix <- History[[iter-1]]$Theta_a
      Sigma_array <- History[[iter-1]]$Sigma_a
    }
    
    ### update U_n ###
    U_n <- rgamma(1, shape=data_size, rate=sum(S_vector))
    
    ### update c ###
    M_old <- History[[iter-1]]$M
    prob_mat <- sapply(1:M_old, function(x) S_vector[x] * dmvnorm(t(data), Theta_matrix[,x], as.matrix(Sigma_array[,,x])))
    c <- apply(prob_mat, 1, function(x) sample(1:M_old, 1, prob = x, replace = TRUE))
    c_table <- tabulate(c, M_old)
    k <- sum(c_table != 0)
    
    ### rename c ###
    ind <- c(which(c_table > 0), which(c_table == 0))
    S_a <- S_vector[ind, drop=FALSE]
    Theta_a <- Theta_matrix[, ind, drop=FALSE]
    Sigma_a <- Sigma_array[, , ind, drop=FALSE]
    
    c_rename <- rep(NA_integer_, data_size)
    for (i in 1:length(ind)) {
      c_rename[c == i] <- which(ind == i)
    }
    c_table <- tabulate(c_rename, k)
    
    ### update M and M_na ###
    M <- sample_M_TS_Ga_Poi(M_max, k, U_n, Ga, lambda, c_rename)
    M_na <- M - k
    
    MH_result_Ga <- sample_Ga_telescope(Ga, df1, df2, U_n, k, M, c_rename, s_pro)
    Ga <- MH_result_Ga$Ga
    acc <- MH_result_Ga$acc
    gamma <- Ga / M
    
    ### update allocated parameters:S_a, Theta_a and Sigma_a ###
    c_table <- table(c_rename)
    for (i in 1:k) {
      S_a[i] <- rgamma(1, shape = c_table[i]+gamma, rate=U_n+1)
    }
    
    c_i <- c_0 + c_table / 2
    Sigma_a_inv <- array(0, dim = c(dimension, dimension, k))
    for (i in 1:k) {
      diff_i <- sweep(data[, which(c_rename == i), drop=FALSE], 1, Theta_a[, i], FUN = "-")
      S_i <- C_0_inv + 0.5 * crossprod(t(diff_i))
      sig_bayesm <- bayesm::rwishart(2 * c_i[i], 0.5 * chol2inv(chol(S_i)))
      Sigma_a[,,i] <- sig_bayesm$IW
      Sigma_a_inv[,,i] <- sig_bayesm$W
    }
    
    for (i in 1:k) {
      N_i <- c_table[i]
      A_i <- chol2inv(chol(B_0_inv + N_i * Sigma_a_inv[,,i]))
      m_i <- A_i %*% ( B_0_inv %*% b_0 + N_i * Sigma_a_inv[,,i] %*% apply(data[, which(c_rename == i), drop=FALSE], 1, mean))
      Theta_a[,i] <- t(chol(as.matrix(A_i))) %*% rnorm(dimension) + m_i
    }
    
    g_K <- g_0 + k * c_0
    C_0_bayesm <- bayesm::rwishart(2*g_K, 0.5*chol2inv(chol(G_0 + rowSums(Sigma_a_inv[, , 1:k, drop = FALSE], dims = 2))))
    C_0 <- C_0_bayesm$IW
    C_0_inv <- C_0_bayesm$W
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na <- rgamma(M_na, shape = gamma, rate = U_n+1)
      Theta_na <- matrix(0, dimension, M_na)
      Sigma_na <- array(0, dim = c(dimension, dimension, M_na))
      for (i in sample(1:M_na)) {
        Sigma_na[,,i] <- bayesm::rwishart(2*c_0, 0.5*chol2inv(chol(C_0)))$IW
        Theta_na[,i] <- rmvnorm(1, b_0, B_0)
      }
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a, S_na=S_na,
                              Theta_a=Theta_a, Theta_na=Theta_na,
                              Sigma_a=Sigma_a, Sigma_na=Sigma_na,
                              C_0 = C_0,
                              lambda=lambda,
                              shape = gamma)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Theta_a=Theta_a,
                              Sigma_a=Sigma_a,
                              C_0 = C_0,
                              lambda=lambda,
                              shape = gamma)
    }
    # cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

TS_Ga_MultiNormal_DMFM_NB <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, r, p, k_ini, M_max){
  ## Model: M-1 \sim NB(r, p), r > 0, p \in (0, 1)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma, M \sim IGau(gamma/M, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
  ##        
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        b_0 : prior mean vector of N_r(b_0, B_0)
  ##        B_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        c_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        g_0, G_0 : hyperparameter for the Wishart prior of C.
  ##        r, p : hyperparameter for the negative binomial prior of M-1
  ##        k_ini : initial number of occupied clusters.
  ##        M_max : upper bound for sampling M
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means, r × k and r × M_na
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        C_0 : hierarchical positive definite symmetric matrix
  ##        shape : shape parameter of h
  History <- vector("list", MCMC_iteration)
  data_size <- length(data[1,])
  dimension <- length(data[,1])
  B_0_inv <- chol2inv(chol(B_0))
  
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
  Ga <- 1.0
  gamma <- Ga
  k <- k_ini
  M_na <- 0
  M <- k + M_na
  c_kmeans <- kmeans(t(data), centers = k, nstart = 100)
  c <- c_kmeans$cluster
  
  C_0 <- g_0 * chol2inv(chol(G_0))
  C_0_inv <- chol2inv(chol(C_0))
  
  S_a <- rep(1/M, k)
  Theta_a <- t(c_kmeans$centers)
  Sigma_a <- array(0, dim = c(dimension, dimension, k))
  Sigma_a[, , 1:k] <- 0.5 * C_0
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, 
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       C_0 = C_0,
                       shape = gamma)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector <- append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
      Theta_matrix <- abind(History[[iter-1]]$Theta_a, History[[iter-1]]$Theta_na)
      Sigma_array <- abind(History[[iter-1]]$Sigma_a, History[[iter-1]]$Sigma_na)
    }else{
      S_vector <- History[[iter-1]]$S_a
      Theta_matrix <- History[[iter-1]]$Theta_a
      Sigma_array <- History[[iter-1]]$Sigma_a
    }
    
    ### update U_n ###
    U_n <- rgamma(1, shape=data_size, rate=sum(S_vector))
    
    ### update c ###
    M_old <- History[[iter-1]]$M
    prob_mat <- sapply(1:M_old, function(x) S_vector[x] * dmvnorm(t(data), Theta_matrix[,x], as.matrix(Sigma_array[,,x])))
    c <- apply(prob_mat, 1, function(x) sample(1:M_old, 1, prob = x, replace = TRUE))
    c_table <- tabulate(c, M_old)
    k <- sum(c_table != 0)
    
    ### rename c ###
    ind <- c(which(c_table > 0), which(c_table == 0))
    S_a <- S_vector[ind, drop=FALSE]
    Theta_a <- Theta_matrix[, ind, drop=FALSE]
    Sigma_a <- Sigma_array[, , ind, drop=FALSE]
    
    c_rename <- rep(NA_integer_, data_size)
    for (i in 1:length(ind)) {
      c_rename[c == i] <- which(ind == i)
    }
    c_table <- tabulate(c_rename, k)
    
    ### update M_na ###
    psi_u <- exp(log_psi_ga(U_n, gamma))
    M <- sample_M_TS_Ga_NB(M_max, k, U_n, Ga, r, p, c_rename)
    M_na <- M - k
    
    ### update gamma ###
    MH_result_Ga <- sample_Ga_telescope(Ga, df1, df2, U_n, k, M, c_rename, s_pro)
    Ga <- MH_result_Ga$Ga
    acc <- MH_result_Ga$acc
    gamma <- Ga / M
    
    ### update allocated parameters:S_a, Theta_a and Sigma_a ###
    for (i in 1:k) {
      S_a[i] <- rgamma(1, shape = c_table[i]+gamma, rate=U_n+1)
    }
    
    c_i <- c_0 + c_table / 2
    Sigma_a_inv <- array(0, dim = c(dimension, dimension, k))
    for (i in 1:k) {
      diff_i <- sweep(data[, which(c_rename == i), drop=FALSE], 1, Theta_a[, i], FUN = "-")
      S_i <- C_0_inv + 0.5 * crossprod(t(diff_i))
      sig_bayesm <- bayesm::rwishart(2 * c_i[i], 0.5 * chol2inv(chol(S_i)))
      Sigma_a[,,i] <- sig_bayesm$IW
      Sigma_a_inv[,,i] <- sig_bayesm$W
    }
    
    for (i in 1:k) {
      N_i <- c_table[i]
      A_i <- chol2inv(chol(B_0_inv + N_i * Sigma_a_inv[,,i]))
      m_i <- A_i %*% ( B_0_inv %*% b_0 + N_i * Sigma_a_inv[,,i] %*% apply(data[, which(c_rename == i), drop=FALSE], 1, mean))
      Theta_a[,i] <- t(chol(as.matrix(A_i))) %*% rnorm(dimension) + m_i
    }
    
    g_K <- g_0 + k * c_0
    C_0_bayesm <- bayesm::rwishart(2*g_K, 0.5*chol2inv(chol(G_0 + rowSums(Sigma_a_inv[, , 1:k, drop = FALSE], dims = 2))))
    C_0 <- C_0_bayesm$IW
    C_0_inv <- C_0_bayesm$W
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na <- rgamma(M_na, shape = gamma, rate = U_n+1)
      Theta_na <- matrix(0, dimension, M_na)
      Sigma_na <- array(0, dim = c(dimension, dimension, M_na))
      for (i in sample(1:M_na)) {
        Sigma_na[,,i] <- bayesm::rwishart(2*c_0, 0.5*chol2inv(chol(C_0)))$IW
        Theta_na[,i] <- rmvnorm(1, b_0, B_0)
      }
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a, S_na=S_na,
                              Theta_a=Theta_a, Theta_na=Theta_na,
                              Sigma_a=Sigma_a, Sigma_na=Sigma_na,
                              C_0 = C_0,
                              shape = gamma)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Theta_a=Theta_a,
                              Sigma_a=Sigma_a,
                              C_0 = C_0,
                              shape = gamma)
    }
    #cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

TS_Ga_MultiNormal_DMFM_BNB <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, al_lam, a_pi, b_pi, k_ini, M_max){
  ## Model: M-1 \sim BNB(al_lam, a_pi, b_pi)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma, M \sim Gamma(gamma/M, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        b_0 : prior mean vector of N_r(b_0, B_0)
  ##        B_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        c_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        g_0, G_0 : hyperparameter for the Wishart prior of C.
  ##        al_lam, a_pi, b_pi : hyperparameter for the beta negative binomial prior of M-1
  ##        k_ini : initial number of occupied clusters.
  ##        Alpha : shape parameter of IGau(Alpha, 1)
  ##        M_max : upper bound for sampling M
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means, r × k and r × M_na
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        C_0 : hierarchical positive definite symmetric matrix
  History <- vector("list", MCMC_iteration)
  data_size <- length(data[1,])
  dimension <- length(data[,1])
  B_0_inv <- chol2inv(chol(B_0))
  
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
  k <- k_ini
  M_na <- 0
  M <- k + M_na
  Ga <- 1.0
  gamma <- Ga
  c_kmeans <- kmeans(t(data), centers = k, nstart = 100)
  c <- c_kmeans$cluster
  
  C_0 <- g_0 * chol2inv(chol(G_0))
  C_0_inv <- chol2inv(chol(C_0))
  
  S_a <- rep(1/M, k)
  Theta_a <- t(c_kmeans$centers)
  Sigma_a <- array(0, dim = c(dimension, dimension, k))
  Sigma_a[, , 1:k] <- 0.5 * C_0
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, 
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       C_0 = C_0,
                       shape = gamma)
  
  ### The MCMC starts. ###
  for (iter in 2:MCMC_iteration) {
    if(History[[iter-1]]$M_na > 0){
      S_vector <- append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
      Theta_matrix <- abind(History[[iter-1]]$Theta_a, History[[iter-1]]$Theta_na)
      Sigma_array <- abind(History[[iter-1]]$Sigma_a, History[[iter-1]]$Sigma_na)
    }else{
      S_vector <- History[[iter-1]]$S_a
      Theta_matrix <- History[[iter-1]]$Theta_a
      Sigma_array <- History[[iter-1]]$Sigma_a
    }
    
    ### update U_n ###
    U_n <- rgamma(1, shape=data_size, rate=sum(S_vector))
    
    ### update c ###
    M_old <- History[[iter-1]]$M
    prob_mat <- sapply(1:M_old, function(x) S_vector[x] * dmvnorm(t(data), Theta_matrix[,x], as.matrix(Sigma_array[,,x])))
    c <- apply(prob_mat, 1, function(x) sample(1:M_old, 1, prob = x, replace = TRUE))
    c_table <- tabulate(c, M_old)
    k <- sum(c_table != 0)
    
    ### rename c ###
    ind <- c(which(c_table > 0), which(c_table == 0))
    S_a <- S_vector[ind, drop=FALSE]
    Theta_a <- Theta_matrix[, ind, drop=FALSE]
    Sigma_a <- Sigma_array[, , ind, drop=FALSE]
    
    c_rename <- rep(NA_integer_, data_size)
    for (i in 1:length(ind)) {
      c_rename[c == i] <- which(ind == i)
    }
    c_table <- tabulate(c_rename, k)
    
    ### update M and M_na ###
    M <- sample_M_TS_Ga_BNB(M_max, k, U_n, Ga, al_lam, a_pi, b_pi, c_rename)
    M_na <- M - k
    
    MH_result_Ga <- sample_Ga_telescope(Ga, df1, df2, U_n, k, M, c_rename, s_pro)
    Ga <- MH_result_Ga$Ga
    acc <- MH_result_Ga$acc
    gamma <- Ga / M
    
    ### update allocated parameters:S_a, Theta_a and Sigma_a ###
    c_table <- table(c_rename)
    for (i in 1:k) {
      S_a[i] <- rgamma(1, shape = c_table[i]+gamma, rate=U_n+1)
    }
    
    c_i <- c_0 + c_table / 2
    Sigma_a_inv <- array(0, dim = c(dimension, dimension, k))
    for (i in 1:k) {
      diff_i <- sweep(data[, which(c_rename == i), drop=FALSE], 1, Theta_a[, i], FUN = "-")
      S_i <- C_0_inv + 0.5 * crossprod(t(diff_i))
      sig_bayesm <- bayesm::rwishart(2 * c_i[i], 0.5 * chol2inv(chol(S_i)))
      Sigma_a[,,i] <- sig_bayesm$IW
      Sigma_a_inv[,,i] <- sig_bayesm$W
    }
    
    for (i in 1:k) {
      N_i <- c_table[i]
      A_i <- chol2inv(chol(B_0_inv + N_i * Sigma_a_inv[,,i]))
      m_i <- A_i %*% ( B_0_inv %*% b_0 + N_i * Sigma_a_inv[,,i] %*% apply(data[, which(c_rename == i), drop=FALSE], 1, mean))
      Theta_a[,i] <- t(chol(as.matrix(A_i))) %*% rnorm(dimension) + m_i
    }
    
    g_K <- g_0 + k * c_0
    C_0_bayesm <- bayesm::rwishart(2*g_K, 0.5*chol2inv(chol(G_0 + rowSums(Sigma_a_inv[, , 1:k, drop = FALSE], dims = 2))))
    C_0 <- C_0_bayesm$IW
    C_0_inv <- C_0_bayesm$W
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na <- rgamma(M_na, shape = gamma, rate = U_n+1)
      Theta_na <- matrix(0, dimension, M_na)
      Sigma_na <- array(0, dim = c(dimension, dimension, M_na))
      for (i in sample(1:M_na)) {
        Sigma_na[,,i] <- bayesm::rwishart(2*c_0, 0.5*chol2inv(chol(C_0)))$IW
        Theta_na[,i] <- rmvnorm(1, b_0, B_0)
      }
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a, S_na=S_na,
                              Theta_a=Theta_a, Theta_na=Theta_na,
                              Sigma_a=Sigma_a, Sigma_na=Sigma_na,
                              C_0 = C_0,
                              shape = gamma)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Theta_a=Theta_a,
                              Sigma_a=Sigma_a,
                              C_0 = C_0,
                              shape = gamma)
    }
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}
