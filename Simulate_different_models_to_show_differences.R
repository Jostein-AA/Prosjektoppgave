library(ggplot2)
library(tidyverse)
require(mgcv)
library(MASS)
library(ar.matrix)
library(INLA)
library(LaplacesDemon)
#----------------------------------------------

#define necessary preliminaries

#Number of samples
n = 10000
T = 1000

#Simulate AR(1)

a = 0.9
tau = 2.5
#for(i in x_axis){
#}

sim_AR1 <- function(x){
  a = 0.9
  tau = 2.5
  T = 1000
  set.seed(5469333)
  increments <- rnorm(n = T, mean = 0, sd = sqrt(1/tau))
  sampled <- rep(0, T)
  sampled[1] = rnorm(n = 1, mean = 0, sd = sqrt(1/(tau * (1 - a**2))))
  for(i in 2:T){
    sampled[i] <- a * sampled[i - 1] +  increments[i]
  }
  return(sampled)
}

#test <- sim_AR1(a, tau, T)
#plot(1:T, test)
ptm <- Sys.time()
ar1_samples <- lapply(1:n, FUN = sim_AR1)
sim_time <- Sys.time() - ptm
print(sim_time)

#Find quantiles for plotting
#first [[i]] indicates which sample, while second [j] indicates time point in that sample 
print(ar1_samples[[1]][1])

#Should maybe reformat?

#Simulate RW(1)




#Plot AR(1) against RW(1)


