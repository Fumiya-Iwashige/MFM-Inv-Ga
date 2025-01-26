rm(list=ls())

#library(extraDistr)
library(GIGrvg)
library(statmod)
library(ggplot2)
library(abind)
library(fossil)

MCMC_iteration = 2000
burn_in = 1000
Alpha = gamma = 1

df_IG = data.frame(
  posterior_mean_of_M = 0,
  error_of_M = 0,
  posteriot_prob_of_M_na = 0,
  posterior_randindex = 0,
  map_randindex = 0,
  posterior_mean_of_lambda = 0,
  posterior_mean_of_U_n = 0
)

pp_M_IG = data.frame(
  "1" = 0,
  "2" = 0,
  "3" = 0,
  "4" = 0,
  "5" = 0,
  "6" = 0,
  "7" = 0,
  "8" = 0,
  "9" = 0,
  "10" = 0
)

pp_k_IG = pp_M_IG

df_Ga = df_IG
pp_M_Ga = pp_M_IG
pp_k_Ga = pp_k_IG

community_number = 3
true_M = community_number
node_size = 150

elements_Q = matrix(0.1, community_number, community_number)
true_Q = matrix(elements_Q, true_M, true_M)
diag(true_Q) = 0.8

set.seed(11111)
data = matrix(0, node_size, node_size)
labels_true = rep(0, node_size)
for (i in 1:node_size) {
  labels_true[i] = i %% true_M
}
labels_true[labels_true == 0] = true_M
labels_true = sort(labels_true)
for (jj in 1:node_size) {
  for (kk in jj:node_size) {
    data[jj, kk] = rbinom(1, 1, prob = true_Q[labels_true[jj], labels_true[kk]])
    data[kk, jj] = data[jj, kk]
  }
}
diag(data) = 0

gamma = 1.0
Alpha = 1.0
a_lambda = 1
b_lambda = 1
Q_shape1 = Q_shape2 = 1
a_lambda = b_lambda = 1

### MFM-Inv-Ga ###
set.seed(11111)
#source("your path including a "IG-functions_network.R" file.")
History_IG = MCMC_IG_network(data,
                             MCMC_iteration,
                             Q_shape1,
                             Q_shape2,
                             a_lambda,
                             b_lambda,
                             Alpha)
df_IG = MCMC_output(MCMC_iteration, burn_in, History_IG, labels_true, true_M)
pp_M_IG = PP_output(MCMC_iteration, burn_in, History_IG, length(pp_M_IG), "M")
pp_k_IG = PP_output(MCMC_iteration, burn_in, History_IG, length(pp_k_IG), "k")


### MFM-Ga ###
set.seed(11111)
#source("your path including a "Ga_functions_network.R" file.")
History_Ga = MCMC_Ga_network(data, 
                             MCMC_iteration, 
                             Q_shape1, 
                             Q_shape2, 
                             a_lambda, 
                             b_lambda, 
                             gamma)
df_Ga = MCMC_output(MCMC_iteration, burn_in, History_Ga, random, true_M)
pp_M_Ga = PP_output(MCMC_iteration, burn_in, History_Ga, length(pp_M_Ga), "M")
pp_k_Ga = PP_output(MCMC_iteration, burn_in, History_Ga, length(pp_k_Ga), "k")
