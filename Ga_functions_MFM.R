BG_Ga_MultiNnormal_MFM_Poi <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, lambda, k_ini){
  ## Model: M-1 \sim Poisson(lambda)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma \sim Gamma(gamma, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
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
  
  gamma <- 1.0
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
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
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

BG_Ga_MultiNnormal_MFM_NB <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, a_lambda, b_lambda, k_ini){
  ## Model: M-1 | lambda \sim Poisson(lambda)
  ##        lambda \sim Gamma(a_lambda, b_lambda)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma \sim IGau(gamma, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
  ##        
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        b_0 : prior mean vector of N_r(b_0, B_0)
  ##        B_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        c_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        g_0, G_0 : hyperparameter for the Wishart prior of C.
  ##        a_lambda, b_lambda : hyperparameter for the Gamma prior of lambda
  ##        k_ini : initial number of occupied clusters.
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means, r × k and r × M_na
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        C_0 : hierarchical positive definite symmetric matrix
  ##        lambda : mean parameter of Poisson distribution
  ##        shape : shape parameter of h
  History <- vector("list", MCMC_iteration)
  data_size <- length(data[1,])
  dimension <- length(data[,1])
  B_0_inv <- chol2inv(chol(B_0))
  
  gamma <- 1.0
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  
  ### initial values ###
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
  
  lambda <- 1.0
  
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
    #cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

BG_Ga_MultiNnormal_MFM_NB_fix <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, a_lambda, b_lambda, gamma, k_ini){
  ## Model: M-1 | lambda \sim Poisson(lambda)
  ##        lambda \sim Gamma(a_lambda, b_lambda)
  ##        y | mu, Sigma \sim N(y | mu, Sigma)
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma \sim IGau(gamma, 1), m=1,...,M
  ##        gamma is fixed.
  ##        
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        b_0 : prior mean vector of N_r(b_0, B_0)
  ##        B_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        c_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        g_0, G_0 : hyperparameter for the Wishart prior of C.
  ##        a_lambda, b_lambda : hyperparameter for the Gamma prior of lambda
  ##        gamma : shape parameter of Gamma(gamma, 1)
  ##        k_ini : initial number of occupied clusters.
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means, r × k and r × M_na
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        C_0 : hierarchical positive definite symmetric matrix
  ##        lambda : mean parameter of Poisson distribution
  History <- vector("list", MCMC_iteration)
  data_size <- length(data[1,])
  dimension <- length(data[,1])
  B_0_inv <- chol2inv(chol(B_0))
  
  ### initial values ###
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
  
  lambda <- 1.0
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, 
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       C_0 = C_0,
                       lambda = lambda)

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
                              lambda=lambda)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Theta_a=Theta_a,
                              Sigma_a=Sigma_a,
                              C_0 = C_0,
                              lambda=lambda)
    }
    #cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

BG_Ga_MultiNnormal_MFM_BNB <- function(data, MCMC_iteration, b_0, B_0, c_0, g_0, G_0, a_lambda, a_pi, b_pi, k_ini){
  ## Model: M-1 | lambda \sim Poisson(lambda)
  ##        lambda | b_lambda \sim Gamma(a_lambda, b_lambda)
  ##        b_lambda \sim Betaprime(a_pi, b_pi)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu \sim N_r(b_0, B_0), Sigma^{-1} | C \sim W(c_0, C), C \sim W(g_0, G_0)
  ##        h_m | gamma \sim Gamma(gamma, 1), m=1,...,M
  ##        gamma \sim F(df1, df2)
  ##        
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        b_0 : prior mean vector of N_r(b_0, B_0)
  ##        B_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        c_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        g_0, G_0 : hyperparameter for the Wishart prior of C.
  ##        a_lambda : hyperparameter for the Gamma prior of lambda
  ##        a_pi, b_pi : hyperparameter for the Betaprime prior of b_lambda
  ##        k_ini : initial number of occupied clusters.
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means, r × k and r × M_na
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        C_0 : hierarchical positive definite symmetric matrix
  ##        lambda : mean parameter of Poisson distribution
  
  History <- vector("list", MCMC_iteration)
  data_size <- length(data[1,])
  dimension <- length(data[,1])
  B_0_inv <- chol2inv(chol(B_0))
  
  df1 <- 1
  df2 <- 1
  s_pro <- 1.0
  acc_blam <- 0
  
  ### initial values ###
  gamma <- 1.0
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
  
  lambda <- 1.0
  b_lambda <- 1.0
  
  History[[1]] <- list(k=k, M_na=M_na, M=M,
                       c=c, S_a=S_a, 
                       Theta_a=Theta_a, Sigma_a=Sigma_a,
                       C_0 = C_0)
  
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
    
    ### update b_lambda ###
    log_b_lambda_pro <- log(b_lambda) + rnorm(1, 0, s_pro)
    b_lambda_pro <- exp(log_b_lambda_pro)
    
    log_accept <- log_post_blambda(b_lambda_pro,a_lambda,lambda,a_pi,b_pi) - log_post_blambda(b_lambda,a_lambda,lambda,a_pi,b_pi) + log(b_lambda_pro) - log(b_lambda)
    accept <- Re(exp(log_accept))
    if (runif(1) <= accept) {
      b_lambda <- b_lambda_pro
      acc_blam <- acc_blam + 1
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
                              C_0 = C_0)
    }else{
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                              c=c_rename, S_a=S_a,
                              Theta_a=Theta_a,
                              Sigma_a=Sigma_a,
                              C_0 = C_0)
    }
    cat("\n", iter, "is over.","\n",table(c_rename))
  }
  
  return(History)
}

BG_Ga_UniNormal_MFM = function(data, MCMC_iteration, m_0, c_0, d_0, D_0, w_0, W_0, a_lambda, b_lambda, gamma){
  ## Model: M-1 | lambda \sim Poisson(lambda)
  ##        lambda \sim Gamma(a_lambda, b_lambda)
  ##        y | mu, Sigma \sim N(y | mu, Sigma) 
  ##        mu | Sigma, tau \sim N_r(m_0, Sigma / B_0), 
  ##        B_0 \sim Gamma(w_0, W_0)
  ##        Sigma^{-1} | C_0 \sim  Gamma(c_0, C_0), 
  ##        C_0 \sim Gamma(d_0, D_0),
  ##        h_m is Gamma(gamma, 1), m=1,...,M
  ##        gamma is fixed
  
  ################################################################
  
  ## Input: data : A numeric matrix of size (r × n). Each column is an observation
  ##        MCMC_iteration : Number of MCMC iterations
  ##        m_0 : prior mean vector of N_r(b_0, B_0)
  ##        c_0 : prior covariance matrix of N_r(b_0, B_0)
  ##        d_0, D_0 : hyperparameter for the Wishart prior of Sigma^{-1} 
  ##        w_0, W_0 : hyperparameter for the Wishart prior of C.
  ##        a_lambda : hyperparameter for the Gamma prior of lambda
  ##        gamma : shape parameter of Gamma(gamma, 1)
  
  ################################################################
  
  ## Return: A list (length MCMC_iteration) containing all parameters:
  ##        M_na : number of empty components
  ##        M : number of components
  ##        c : label vector
  ##        S_a, S_na: un-normalized component weights
  ##        Theta_a, Theta_na : component means
  ##        Sigma_a, Sigma_na : component covariance matrices
  ##        lambda : mean parameter of Poisson distribution
  ##        B_0 : smoothing parameter
  ##        C_0 : rate parameter of Gamma distribution for Sigma
  
  History <- vector("list", MCMC_iteration)
  data_size <- length(data)
  
  ### initial values ###
  K <- 10
  k <- sample(1:K, 1, replace = TRUE, prob = rep(1/K, K))
  M_na <- sample(1:K, 1, replace = TRUE, prob = rep(1/K, K))
  M <- k + M_na
  c <- sample(1:k, data_size, replace = TRUE, prob = rep(1/k, k)); c[1:k] = 1:k 
  S_a <- rep(1, k)
  Theta_a <- rnorm(k, 0, 1)
  Sigma_a <- rinvgamma(k, 1, 1)
  S_na <- rep(1, M_na)
  Theta_na <- rnorm(M_na, 0, 1)
  Sigma_na <- rinvgamma(M_na, 0, 1)
  X <- matrix(1, 1, data_size)
  lambda <- rgamma(1, a_lambda, b_lambda)
  B_0 <- rgamma(1, 1, 1)
  C_0 <- rgamma(1, 1, 1)
  
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
      S_vector <- append(History[[iter-1]]$S_a, History[[iter-1]]$S_na)
      Theta_vector <- append(History[[iter-1]]$Theta_a, History[[iter-1]]$Theta_na)
      Sigma_vector <- append(History[[iter-1]]$Sigma_a, History[[iter-1]]$Sigma_na)
    }else{
      S_vector <- History[[iter-1]]$S_a
      Theta_vector <- History[[iter-1]]$Theta_a
      Sigma_vector <- History[[iter-1]]$Sigma_a
    }
    
    ### update U_n ###
    U_n <- rgamma(1, shape=data_size, rate=sum(S_vector))
    
    ### update c ###
    c_old <- History[[iter-1]]$c
    M_old <- History[[iter-1]]$M
    c_update <- rep(0, data_size)
    for (i in sample(1:data_size)) {
      weight_i <- rep(0, M_old)
      for (j in 1:M_old) {
        weight_i[j] <- S_vector[j]*dnorm(data[i], Theta_vector[j], sqrt(Sigma_vector[j]), log=FALSE)
      }
      if(sum(weight_i) == 0){
        for (j in 1:M_old) {
          weight_i[j] <- log(S_vector[j]) + dnorm(data[i], Theta_vector[j], sqrt(Sigma_vector[j]), log=TRUE)
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
    Theta_a <- Theta_vector[c_unique]
    Sigma_a <- Sigma_vector[c_unique]
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
    
    ### update allocated parameters:S_a, Theta_a and Sigma_a ###
    c_table <- table(c_rename)
    for (i in sample(1:k)) {
      S_a[i] <- rgamma(1, shape = c_table[i]+gamma, rate=U_n+1)
      
      N_i <- c_table[i]
      Y_i <- data[which(c_rename == i)]
      X_i <- X[which(c_rename == i)]
      B_i <- 1 / (N_i + B_0)
      M_i <- (sum(Y_i) + B_0 * m_0) * B_i
      c_i <- c_0 + N_i/2
      C_k <- B_0 * (Theta_a[i] - m_0)^2 + sum((Y_i - M_i)^2)
      
      Theta_a[i] <- rnorm(1, M_i, sqrt(Sigma_a[i]*B_i))
      Sigma_a[i] <- rinvgamma(1, c_i, C_0 + C_k / 2)
    }
    
    ### update unallocated parameters and put result###
    if(M_na > 0){
      S_na <- rgamma(M_na, shape = gamma, rate = U_n+1)
      Theta_na <- rep(0, M_na)
      Sigma_na <- rep(0, M_na)
      sample_vector <- sample(1:M_na)
      for (i in sample_vector) {
        Sigma_na[i] <- rinvgamma(1, c_0, C_0)
        Theta_na[i] <- rnorm(1, m_0, sqrt(Sigma_na[i]/B_0))
      }
      
      ### update B_0 and C_0###
      W_1 <- (k+M_na) /2 + w_0
      W_2 <- W_0 + sum((append(Theta_a, Theta_na) - m_0)^2 / append(Sigma_a, Sigma_na)) * 0.5
      B_0 <- rgamma(1, W_1, W_2)
      
      D_1 <- (k + M_na) * c_0 + d_0
      D_2 <- sum(1/append(Sigma_a, Sigma_na)) + D_0
      C_0 <- rgamma(1, D_1, D_2)
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
                             c=c_rename, S_a=S_a, S_na=S_na,
                             Theta_a=Theta_a, Theta_na=Theta_na,
                             Sigma_a=Sigma_a, Sigma_na=Sigma_na,
                             lambda=lambda,
                             B_0 = B_0,
                             C_0 = C_0)
    }else{
      ### update B_0 ###
      W_1 <- k / 2 + w_0
      W_2 <- W_0 + sum((Theta_a - m_0)^2 / Sigma_a) * 0.5
      B_0 <- rgamma(1, W_1, W_2)
      
      D_1 <- (k + M_na) * c_0 + d_0
      D_2 <- sum(1/Sigma_a) + D_0
      C_0 <- rgamma(1, D_1, D_2)
      
      History[[iter]] <- list(k=k, M_na=M_na, M=k+M_na, U = U_n,
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
