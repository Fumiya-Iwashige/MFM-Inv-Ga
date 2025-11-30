rm(list=ls())

source("IG_functions_MFM.R")
source("IG_functions_DMFM.R")
source("IG_functions_general.R")
source("Ga_functions_MFM.R")
source("Ga_functions_DMFM.R")
source("Ga_functions_general.R")
source("summary.R")

library(abind)
library(bmixture)
library(ggplot2)
library(LaplacesDemon)
library(statmod)
library(GIGrvg)
library(reshape2)

repetition <- 5
MCMC_iteration <- 10000
burn_in <- 9000

density_IG <- vector("list", repetition)

### data setting ###
data("galaxy")
data <- galaxy
data_size <- length(data)
data_length <- max(data) - min(data)
grid_y <- seq(min(data),max(data), len=1000)

### hyper parameters ###
c_0 <- 2; C_0 <- 1
w_0 <- 1/2; W_0 <- 100/2
d_0 <- 0.2; D_0 <- 10 / (data_length^2)
a_lambda <- 1.0
b_lambda <- 1 / 5
m_0 <- data_length / 2


### MFM-IGau ###
set.seed(12345)
Alpha <- 1.0
History_IG <- BG_IG_UniNormal_MFM(data,
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

### density estimate ###
estimate_IG <- density_estimate(grid_y, History_IG, MCMC_iteration, burn_in)
g <- ggplot(estimate_IG, aes(x = x, y = posterior_mean))
g <- g + geom_histogram(data = data.frame(galaxy), aes(x = galaxy, y = ..density..), bins=30, fill="white", colour="black") +
  geom_ribbon(aes(ymin = quontile_0025, ymax = quontile_0975), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = posterior_mean), color = "blue", size = 0.5) + 
  xlab("velocity") + 
  ylab("density") + 
  ggtitle("density estimation via MFM-IGau") +
  theme_light()
plot(g)


### MFM-Ga ###
set.seed(12345)
gamma <- 1.0
History_Ga <- BG_Ga_UniNormal_MFM(data, 
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

### density estimate ###
estimate_Ga <- density_estimate(grid_y, History_Ga, MCMC_iteration, burn_in)

g <- ggplot(estimate_Ga, aes(x = x, y = posterior_mean))
g <- g + geom_histogram(data = data.frame(galaxy), aes(x = galaxy, y = ..density..), bins=30, fill="white", colour="black") +
  geom_ribbon(aes(ymin = quontile_0025, ymax = quontile_0975), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = posterior_mean), color = "blue", size = 0.5) + 
  xlab("velocity") + 
  ylab("density") + 
  ggtitle("density estimation via MFM-Ga") +
  theme_light()
plot(g)