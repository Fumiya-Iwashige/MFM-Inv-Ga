library(abind)
library(dplyr)
library(LaplacesDemon)
library(statmod)
library(mvtnorm)
library(GIGrvg)
library(extraDistr)
library(fossil)
library(tidyr)

rm(list=ls())
source("IG_functions_MFM.R")
source("IG_functions_DMFM.R")
source("IG_functions_general.R")
source("Ga_functions_MFM.R")
source("Ga_functions_DMFM.R")
source("Ga_functions_general.R")
source("summary.R")

data_size <- 150
MCMC_iteration <- 30000
burn_in <- 20000

cluster_number <- 4
true_M <- cluster_number
dimension <- 2
components_ratio <- rep(1/true_M, true_M)

### mu is designed the same as in Schnatter et al. (2021)
true_mu <- matrix(0, dimension, cluster_number)
true_1 <- c(2, 6, 10, 14); true_2 <- c(0, 5)
true_mu[1,1] <- true_1[1]; true_mu[2,1] <- true_2[1]
true_mu[1,2] <- true_1[1]; true_mu[2,2] <- true_2[2]
true_mu[1,3] <- true_1[2]; true_mu[2,3] <- true_2[1]
true_mu[1,4] <- true_1[2]; true_mu[2,4] <- true_2[2]

sig <- 0.5
true_sigma <- array(0, dim = c(dimension, dimension, cluster_number))
for (i in 1:cluster_number) {
  true_sigma[,,i] <- diag(1.0, dimension) * sig
}

true_labels <- sample(1:cluster_number, data_size, replace = TRUE, prob = components_ratio)
data_set <- matrix(NA_real_, nrow = dimension, ncol = data_size)
for (i in 1:data_size) {
  data_set[,i] <- rmvnorm(1, mean = true_mu[,true_labels[i]], sigma = true_sigma[,,true_labels[i]])
}

### hyper parameters independent on data ###
c_0 <- dimension + 1.5
a_lambda <- 1.0
b_lambda <- 1.0
k_ini <- 10

R <- apply(t(data_set), 2, function(x) diff(range(x)))
b_0 <- apply(t(data_set), 2, median)
B0 <- rep(1, dimension)
B_0 <- diag((R^2) * B0)
c_0 <- 2.5 + (dimension-1)/2
g_0 <- 0.5 + (dimension-1)/2
G_0 <- 100 * g_0/c_0 * diag((1/R^2), nrow = dimension)

History <- BG_IG_MultiNnormal_MFM_BNB(data_set,
                                      MCMC_iteration,
                                      b_0,
                                      B_0,
                                      c_0,
                                      g_0,
                                      G_0,
                                      a_lambda,
                                      a_pi=4,
                                      b_pi=3,
                                      k_ini)
MCMC_output(MCMC_iteration, burn_in, History, true_labels, "MFM-IG-BNB")

History <- TS_IG_MultiNormal_DMFM_BNB(data_set,
                                      MCMC_iteration,
                                      b_0,
                                      B_0,
                                      c_0,
                                      g_0,
                                      G_0,
                                      al_lam=1,
                                      a_pi=4,
                                      b_pi=3,
                                      k_ini,
                                      M_max=100)
MCMC_output(MCMC_iteration, burn_in, History, true_labels, "DMFM-IG-BNB")

History <- BG_Ga_MultiNnormal_MFM_BNB(data_set, 
                                      MCMC_iteration, 
                                      b_0, 
                                      B_0, 
                                      c_0, 
                                      g_0, 
                                      G_0, 
                                      a_lambda, 
                                      a_pi=4, 
                                      b_pi=3, 
                                      k_ini)
MCMC_output(MCMC_iteration, burn_in, History, true_labels, "MFM-Ga-BNB")

History <- TS_Ga_MultiNormal_DMFM_BNB(data_set,
                                     MCMC_iteration,
                                     b_0,
                                     B_0,
                                     c_0,
                                     g_0,
                                     G_0,
                                     al_lam=1, 
                                     a_pi=4, 
                                     b_pi=3,
                                     k_ini,
                                     M_max=100)
MCMC_output(MCMC_iteration, burn_in, History, true_labels, "DMFM-Ga-BNB")
