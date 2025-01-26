rm(list=ls())

library(abind)
library(bmixture)
library(ggplot2)
library(LaplacesDemon)
library(statmod)
library(GIGrvg)

MCMC_iteration = 20000
burn_in = 10000

df_IG = data.frame(
  posterior_mean_of_M = 0,
  error_of_M = 0,
  posteriot_prob_of_M_na = 0,
  posterior_randindex = 0,
  map_randindex = 0,
  posterior_mean_of_lambda = 0,
  posterior_mean_of_U_n = 0,
  posterior_mean_of_B_0 = 0,
  posterior_mean_of_C_0 = 0
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
pp_M_Ga = pp_M_IG
pp_k_Ga = pp_k_IG

### data setting ###
data("galaxy")
data = galaxy
data_size = length(data)
grid_y = seq(min(data),max(data), len=1000)

### hyper parameters ###
Alpha = 1.0
gamma = 1.0
data_length = max(data) - min(data)
c_0 = 2; C_0 = 1
w_0 = 1/2; W_0 = 100/2
d_0 = 0.2; D_0 = 10 / (data_length^2)
a_lambda = 1.0
b_lambda = 1 / 5
m_0 = data_length / 2

### MFM-Inv-Ga ###
set.seed(11111)
#source("your path including a "IG_functions.R" file.")
#source("your path including a "density_output.R" file.")
History_IG = MCMC_IG_blocking_uni_normal(data,
                                          MCMC_iteration,
                                          m_0,
                                          c_0,
                                          d_0,
                                          D_0,
                                          w_0,
                                          W_0,
                                          a_lambda,
                                          b_lambda,
                                          Alpha)
df_IG = MCMC_output_density(MCMC_iteration, burn_in, History_IG)
pp_M_IG = PP_output_density(MCMC_iteration, burn_in, History_IG, length(pp_M_IG), "M")
pp_k_IG = PP_output_density(MCMC_iteration, burn_in, History_IG, length(pp_k_IG), "k")
estimate_IG = density_estimate(grid_y, History_IG, MCMC_iteration, burn_in)
density_IG = list(estimate = estimate_IG)

g = ggplot(density_IG$estimate, aes(x = x, y = posterior_mean))
g = g + geom_histogram(data = data.frame(galaxy), aes(x = galaxy, y = ..density..), bins=30, fill="white", colour="black") +
  geom_ribbon(aes(ymin = quontile_0025, ymax = quontile_0975), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = posterior_mean), color = "blue", size = 0.5) + 
  xlab("velocity") + 
  ylab("density") + 
  theme_light()
plot(g)

### MFM-Ga ###
set.seed(11111)
#source("your path including a "Ga_functions.R" file.")
History_Ga = MCMC_Ga_blocking_uni_normal(data, 
                                          MCMC_iteration, 
                                          m_0, 
                                          c_0, 
                                          d_0,
                                          D_0,
                                          w_0, 
                                          W_0, 
                                          a_lambda, 
                                          b_lambda, 
                                          gamma)
df_Ga = MCMC_output_density(MCMC_iteration, burn_in, History_Ga)
pp_M_Ga = PP_output_density(MCMC_iteration, burn_in, History_Ga, length(pp_M_Ga), "M")
pp_k_Ga = PP_output_density(MCMC_iteration, burn_in, History_Ga, length(pp_k_Ga), "k")
estimate_Ga = density_estimate(grid_y, History_Ga, MCMC_iteration, burn_in)
density_Ga = list(estimate = estimate_Ga)

gg = ggplot(density_Ga$estimate, aes(x = x, y = posterior_mean))
gg = gg + geom_histogram(data = data.frame(galaxy), aes(x = galaxy, y = ..density..), bins=30, fill="white", colour="black") +
  geom_ribbon(aes(ymin = quontile_0025, ymax = quontile_0975), fill = "red", alpha = 0.2) +
  geom_line(aes(y = posterior_mean), color = "red", size = 0.5) + 
  xlab("velocity") + 
  ylab("density") + 
  theme_light()
plot(gg)
