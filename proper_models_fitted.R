#Clear environment
rm(list = ls())

#Load necessary libraries
library(INLA)
library(spData)
library(sf)
library(spdep)

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
#This is done for the sake of the space time interaction in this model
ohio_df = ohio_df[order(ohio_df$county, decreasing = F), ] 
ohio_df$space.time <- 1:(n *T); rownames(ohio_df) <- 1:nrow(ohio_df) 

#Make copies of county and year.
#It is done because both year and county is used twice in some formulas
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

#Define Temporal hyperparameters and corresponding priors
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.01)), #Magic numbers
                 rho = list(prior = 'normal', 
                            param = c(0, 0.25)), #Magic numbers
                 mean = list(prior = 'normal',
                             param = c(0, 2.5))) #Magic numbers


#Define Spatial hyperparameters and corresponding priors
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01)), #Magic numbers
                     lambda = list(prior = 'gaussian', 
                                   param = c(0, 1))) #Magic numbers




###
#Define a baseline proper formula: intercept, fixed temporal effect, ar1 temporal and proper-besag spatial (no interactions)
proper_base_formula <- deaths ~ 1 + year +
                                f(year.copy,
                                  model = "ar1",
                                  hyper = ar1_hyper) + 
                                f(county, 
                                  model = "besagproper2",
                                  graph = Besag_prec,
                                  hyper = spatial_hyper)



ptm <- Sys.time()
proper_base_fit <- inla(proper_base_formula, 
                        data = ohio_df,
                        family = "poisson",
                        E = pop_at_risk, 
                        control.compute = list(config = TRUE, # To see constraints later
                                               cpo = T,   # For model selection
                                               waic = T)) # For model selection

time_proper_base = Sys.time()-ptm
print(c("Basic model fitted in: ", time_proper_base))
print(mean(-log(proper_base_fit$cpo$cpo)))

###
#Define a model with intercept, fixed temporal effect and spatiotemporal interaction by ar1 and properbesag
proper_interaction_formula <- deaths ~ 1 + year + 
                                       f(county, 
                                         model = "besagproper2",
                                         graph = Besag_prec,
                                         group = year, 
                                         control.group = list(model = "ar1"))


ptm <- Sys.time()
proper_interaction_fit <- inla(proper_interaction_formula, 
                               data = ohio_df,
                               family = "poisson",
                               E = pop_at_risk, 
                               control.compute = list(config = TRUE, # To see constraints later
                                                      cpo = T,   # For model selection
                                                      waic = T)) # For model selection

time_proper_interaction = Sys.time()-ptm
print(c("Proper Interaction model fitted in: ", time_proper_interaction))
print(mean(-log(proper_interaction_fit$cpo$cpo)))

###
#Define a model with interceot, fixed temporal effect, proper random effects, and proper interaction
proper_full_formula <- deaths ~ 1 + year + 
                                f(year.copy,
                                  model = "ar1",
                                  hyper = ar1_hyper) +
                                f(county, 
                                  model = "besagproper2",
                                  graph = Besag_prec,
                                  hyper = spatial_hyper) + 
                                f(county.copy, 
                                  model = "besagproper2",
                                  graph = Besag_prec,
                                  hyper = spatial_hyper,
                                  group = year, 
                                  control.group = list(model = "ar1"))                                



ptm <- Sys.time()
proper_full_fit <- inla(proper_full_formula, 
                        data = ohio_df,
                        family = "poisson",
                        E = pop_at_risk, 
                        control.compute = list(config = TRUE, # To see constraints later
                                               cpo = T,   # For model selection
                                               waic = T)) # For model selection

time_proper_full = Sys.time()-ptm
print(c("Proper Full model fitted in: ", time_proper_full))
print(mean(-log(proper_full_fit$cpo$cpo)))


#####
#Save INLA objects
ohio_df_changed = ohio_df
save(n, T, ohio_map, ohio_df_changed,
     proper_base_fit, time_proper_base, 
     proper_interaction_fit, time_proper_interaction,
     proper_full_fit, time_proper_full,
     file = "proper_fitted.RData")





