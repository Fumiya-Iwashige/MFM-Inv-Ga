library(GIGrvg)
library(statmod)
library(abind)
library(fossil)
library(igraph)
library(extraDistr)

rm(list=ls())

source("IG_functions_MFM_network.R")
source("IG_functions_DMFM_network.R")
source("IG_functions_general.R")
source("Ga_functions_MFM_network.R")
source("Ga_functions_DMFM_network.R")
source("Ga_functions_general.R")
source("Geng_MFM.R")
source("summary.R")

MCMC_iteration <- 10000
burn_in <- 5000

### data generating ###
community_number <- 2
true_M <- community_number
node_size <- 100

elements_Q <- matrix(0.1, community_number, community_number)
true_Q <- matrix(elements_Q, true_M, true_M)
diag(true_Q) <- 0.8

set.seed(12345)
true_allocation <- rep(0, node_size)
for (i in 1:node_size) {
  true_allocation[i] <- i %% true_M
}
true_allocation[true_allocation == 0] <- true_M

data_set <- matrix(0, nrow = node_size, ncol = node_size)
data_upper <- matrix(0, nrow = node_size, ncol = node_size)
for (i in 1:node_size) {
  for (j in i:node_size) {
    data_set[i, j] <- rbinom(1, 1, prob = true_Q[true_allocation[i], true_allocation[j]])
    data_set[j, i] <- data_set[i, j]
    data_upper[i, j] <- data_set[i, j]
  }
}
diag(data_set) <- 0
diag(data_upper) <- 0

### hypperparameters ###
Q_shape1 <- Q_shape2 <- 1
a_lambda <- b_lambda <- 1

### MFM-IGau ###
History_IG_MFM <- BG_IG_network_MFM_BNB(data_set,
                                 MCMC_iteration,
                                 Q_shape1,
                                 Q_shape2,
                                 a_lambda = 1,
                                 a_pi = 4,
                                 b_pi = 3)
MCMC_output(MCMC_iteration, burn_in, History_IG_MFM, true_allocation, "MFM-IG-BNB")

### DMFM-IGau ###
History_IG_DMFM <- TS_IG_network_DMFM_BNB(data_set, 
                                  MCMC_iteration, 
                                  Q_shape1, 
                                  Q_shape2, 
                                  al_lam=1, 
                                  a_pi=4, 
                                  b_pi=3, 
                                  M_max=node_size)
MCMC_output(MCMC_iteration, burn_in, History_IG_DMFM, true_allocation, "DMFM-IG-NB")

### DMFM-Ga ###
History_Ga_DMFM <- TS_Ga_network_DMFM_BNB(data_set, 
                                  MCMC_iteration, 
                                  Q_shape1, 
                                  Q_shape2, 
                                  al_lam=1, 
                                  a_pi=4, 
                                  b_pi=3, 
                                  M_max=node_size)
MCMC_output(MCMC_iteration, burn_in, History_Ga_DMFM, true_allocation, "DMFM-Ga-NB")