library(ggplot2)
library(tidyverse)
require(mgcv)
library(MASS)
library(ar.matrix)
library(INLA)
library(LaplacesDemon)
#----------------------------------------------

n_samples = 1
dim_ar_rw = 1000

#Simluate AR(1)

#Define different values of a in x_t = ax_(t-1)+epsilon
#a = seq(-0.99, 0.99, length.out = 10)

#Define pairwise observation std deviation
sigma = 1

#Define precision matrix of AR1
#prec_AR1 <- Q.AR1(M = 20, sigma = sigma, rho = 0.7,
#                  sparse = TRUE, vcov =  FALSE)


#Simulate values from AR1
samples <- r.AR1(n = 1, M = dim_ar_rw, sigma = sigma, rho = 0.97)


# move objects to a data frame
ar1_df <- data.frame(obs=c(t(samples)), realization=rep(1:1, each=dim_ar_rw),
                     time=rep(1:dim_ar_rw, 1))

# plot each realization
ggplot(data=ar1_df, aes(time, obs)) +
  geom_line()



#----------------------------------------------
#Simlate AR(2)


AR2.sim <- function(n, ar1, ar2, sd)
{
  Xt = c(0,0)
  for (i in 2:(n-1))
  {
    Xt[i+1] = ar1*Xt[i] + ar2*Xt[i-1] + rnorm(1,mean=0,sd=sd)
  }
  Xt
}

#Poles for z^2-z*a1-a2 = 0. Poles must lie strictly inside unit circle
a1 = 1/2
a2 = 1/4

Xt = AR2.sim(dim_ar_rw, a1, a2, sigma)

ar2_df <- data.frame(obs=c(t(Xt)), realization=rep(1:1, each=dim_ar_rw),
                     time=rep(1:dim_ar_rw, 1))

ggplot(data=ar2_df, aes(time, obs)) +
  geom_line()


#----------------------------------------------
#Simlate RW(1)

#struct_RW1 <- INLA:::inla.rw(n = dim_ar_rw, order = 1, 
#                             scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
#                             sparse = FALSE)


#samples_rw1 = rmvnp(n = 1, mu = rep(0, dim_ar_rw), Omega = struct_RW1)
#samples_rw1 <- r.RW1()


#----------------------------------------------
#Simulate RW(2)

struct_RW2 <- INLA:::inla.rw(n = 1000, order = 2, 
                             scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                             sparse = TRUE)



#----------------------------------------------
#Simulate CAR


#----------------------------------------------
#Simulate ICAR

