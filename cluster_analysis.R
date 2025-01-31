rm(list=ls())
rm(list=ls())

library(abind)
library(LaplacesDemon)
library(statmod)
library(mvtnorm)
library(GIGrvg)
library(fossil)

source("IG_functions.R")
source("Ga_functions.R")

MCMC_iteration = 2000
burn_in = 1000

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

cluster_number = 3
true_M = cluster_number
data_size = 300
dimension = 2
components_ratio = c(0.8, 0.1, 0.1)
#components_ratio = rep(1/true_M, true_M)

true_mu = matrix(0, dimension, cluster_number)
true_1 = c(0, 7.5); true_2 = c(0, 10)
true_mu[1,1] = true_1[1]; true_mu[2,1] = true_2[1]
true_mu[1,2] = true_1[1]; true_mu[2,2] = true_2[2]
true_mu[1,3] = true_1[2]; true_mu[2,3] = true_2[1]
true_mu = true_mu 

true_sigma = diag(1.00, dimension)

set.seed(12345)
data = matrix(0, dimension, data_size)
true_allocation = sample(1:cluster_number, data_size, replace = TRUE, prob = components_ratio)

for (i in 1:data_size) {
  data[, i] = rmvnorm(1, mean = true_mu[,true_allocation[i]], sigma = true_sigma)
}

plot(data[1,], data[2,])

gamma = 1.0
Alpha = 1.0
c_0 = dimension + 1.5
a_lambda = 1
b_lambda = 1
m_0 = apply(data, 1, mean)
C_0 = cov(t(data))

### MFM-Inv-IG ###
History = MCMC_IG_blocking_multi_normal(data, 
                                         MCMC_iteration, 
                                         m_0, 
                                         c_0, 
                                         C_0,
                                         a_lambda, 
                                         b_lambda, 
                                         Alpha)
df_IG = MCMC_output(MCMC_iteration, burn_in, History, true_allocation, true_M)
pp_M_IG = PP_output(MCMC_iteration, burn_in, History, length(pp_M_IG), "M")
pp_k_IG = PP_output(MCMC_iteration, burn_in, History, length(pp_M_IG), "k")

### MFM-Ga ###
History = MCMC_Ga_blocking_multi_normal(data, 
                                        MCMC_iteration, 
                                        m_0, 
                                        c_0, 
                                        C_0,
                                        a_lambda, 
                                        b_lambda, 
                                        gamma)
df_Ga = MCMC_output(MCMC_iteration, burn_in, History, true_allocation, true_M)
pp_M_Ga = PP_output(MCMC_iteration, burn_in, History, length(pp_M_Ga), "M")
pp_k_Ga = PP_output(MCMC_iteration, burn_in, History, length(pp_M_Ga), "k")
