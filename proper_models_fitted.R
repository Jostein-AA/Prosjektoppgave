#Clear environment
rm(list = ls())

#Load necessary libraries
library(INLA)
library(spData)
library(sf)
library(spdep)

#Load utility functions
#source("utilities.R")
#####
#Load data
ohio_df <- read.csv("ohio_df.csv")
ohio_df$year <- ohio_df$year - min(ohio_df$year) + 1

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]

#Get number of counties (n) and number of years (T)
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points

#Firstly, switch ordering from year over county, to county over year
ohio_df = ohio_df[order(ohio_df$county, decreasing = F), ] #This is done for the sake of the space time interaction
ohio_df$space.time <- 1:(n *T); rownames(ohio_df) <- 1:nrow(ohio_df) #in this modell
ohio_df$county.copy <- ohio_df$county
ohio_df$year.copy <- ohio_df$year

#Make precision matrix for spatial Besag
nb <- spdep::poly2nb(ohio_map, queen = FALSE)
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
Besag_prec <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

#####
#Define hyperparameters and corresponding priors
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), #Magic numbers
                 rho = list(prior = 'normal', 
                            param = c(0, 0.25)), #Magic numbers
                 mean = list(prior = 'normal',
                             param = c(0, 2.5))) #Magic numbers


#Spatial hyperparameters 
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01)), #Magic numbers
                     lambda = list(prior = 'gaussian', 
                                   param = c(0, 1))) #Magic numbers

#####
#Define the a basic proper formula 
#(no space-time interaction, but fixed temporal effect and  proper random effects)
proper_formula <- deaths ~ 1 + year + 
                           f(year.copy, 
                             model = "ar1",
                             hyper = ar1_hyper) +
                           f(county, 
                             model = "besagproper2",
                             graph = Besag_prec,
                             group = year, 
                             control.group = list(model = "ar1"))

ptm <- Sys.time()
proper_fit <- inla(proper_formula, 
                   data = ohio_df,
                   family = "poisson",
                   E = pop_at_risk, 
                   control.compute = list(config = TRUE, # To see constraints later
                                          cpo = T,   # For model selection
                                          waic = T)) # For model selection

time_proper = Sys.time()-ptm
print(c("Basic model fitted in: ", time_proper))
print(mean(-log(proper_fit$cpo$cpo)))


#####
#Define proper formula using only fixed temporal effect
# and spatio-temporal interaction
proper_formula_2 <- deaths ~ 1 + year + #year is fixed temporal effect
                              f(county, 
                                model = "besagproper2",
                                graph = Besag_prec,
                                group = year, 
                                control.group = list(model = "ar1"))



ptm <- Sys.time()
proper_fit_2 <- inla(proper_formula_2, 
                   data = ohio_df,
                   family = "poisson",
                   E = pop_at_risk, 
                   control.compute = list(config = TRUE, # To see constraints later
                                          cpo = T,   # For model selection
                                          waic = T)) # For model selection

time_proper_2 = Sys.time()-ptm
print(c("Basic model fitted in: ", time_proper_2))
print(mean(-log(proper_fit_2$cpo$cpo)))


#####
#Define model having fixed temporal effect, 
#random spatial and temporal effects and spatio-temporal interactions all proper
proper_formula_3 <- deaths ~ 1 + year +
                          f(year.copy, #f(year) is temporal random effect (AR(1))
                             model = "ar1",
                             hyper = ar1_hyper) +
                          f(county, 
                            model = "besagproper2", #f(county) is spatial random effect (besagproper)
                            graph = Besag_prec,
                            hyper = spatial_hyper) +
                          f(county.copy,    #f(county.copy, group = year) is spatio-temporal interaction
                            model = "besagproper2",
                            graph = Besag_prec,
                            group = year, 
                            control.group = list(model = "ar1"))
  

ptm <- Sys.time()
proper_fit_3 <- inla(proper_formula_3, 
                     data = ohio_df,
                     family = "poisson",
                     E = pop_at_risk, 
                     control.compute = list(config = TRUE, # To see constraints later
                                            cpo = T,   # For model selection
                                            waic = T)) # For model selection

time_proper_3 = Sys.time()-ptm
print(c("Basic model fitted in: ", time_proper_3))
print(mean(-log(proper_fit_3$cpo$cpo)))

#####
#Save INLA objects
ohio_df_changed = ohio_df
save(n, T, ohio_map, ohio_df_changed,
     proper_fit, time_proper, 
     proper_fit_2, time_proper_2,
     proper_fit_3, time_proper_3,
     file = "proper_fitted.RData")





