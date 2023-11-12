#Clear environment
rm(list = ls())

#Load necessary libraries
library(INLA)
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)
library(ggspatial)
library(gridExtra)
library(ggpubr)


#Set working directory
if(getwd() != "C:/Users/joste/Documents/H2023/Code/Prosjektoppgave"){
  setwd("H2023/Code/Prosjektoppgave/")
}
source("utilities.R")

#Load data
ohio_df <- read.csv("ohio_df.csv")

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]

#Get number of counties (n) and number of years (T)
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points

#Get precision matricies for RW1 and RW2
RW1_prec <- INLA:::inla.rw(n = T, order = 1, 
                             scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                             sparse = TRUE)

RW2_prec <- INLA:::inla.rw(n = T, order = 2, 
                             scale.model = FALSE,
                             sparse = TRUE)

#Extract adjacency structure and create precision matrix for ICAR
nb <- spdep::poly2nb(ohio_map, queen = FALSE)
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
ICAR_prec <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

#Specify hyperparameters with corresponding priors
#Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec.unstruct = list(prior = 'pc.prec', 
                                           param = c(1, 0.01)), #Magic numbers
                      prec.spatial = list(prior = 'pc.prec', 
                                          param = c(1, 0.01)) #Magic numbers
) 

#Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec.unstruct = list(prior = 'pc.prec', 
                                          param = c(1, 0.01)), #Magic numbers
                     prec.spatial = list(prior = 'pc.prec', 
                                         param = c(1, 0.01)) #Magic numbers
)
#Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec",
                                    param=c(1,0.01)))

#Make the base formula
base_formula <- deaths ~ 1 + f(year, 
                               model = 'bym',
                               scale.model = T, 
                               constr = T, 
                               rankdef = 1,
                               graph = RW1_prec,
                               hyper = temporal_hyper) + 
                             f(county, 
                               model = 'bym',
                               scale.model = T,
                               constr = T,
                               rankdef = 1,
                               graph = ICAR_prec,
                               hyper = spatial_hyper)

###

#Fit the base formula 
ptm <- Sys.time()
basic_model_fit <- inla(base_formula, data = ohio_df, family = "poisson",
                        E = pop_at_risk, 
                        control.compute = list(config = TRUE, # To see constraints later
                                               cpo = T,   # For model selection
                                               waic = T), # For model selection
                        ) 
time_base = Sys.time()-ptm

#Inference on the basic model
print_cpo_etc(basic_model_fit, time_base)

#See a lot of things (intercept, precisions, effects, etc)
plot(basic_model_fit)

#Plot the intercept
plot_intercept(basic_model_fit)

#Plot posterior distributions of precision of random effects
plot_precisions_random_effects(basic_model_fit)

#See Structured temporal effect
plot_temporal_effect(basic_model_fit)

#See structured spatial effect plotted as heatmap
plot_spatial_effect(ohio_map, basic_model_fit)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(basic_model_fit)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, basic_model_fit)

###

#Update base formula to also contain iid interaction
typeI_formula <- update(base_formula,  ~. + f(space_time_unstructured,
                                              model="iid", #Has to be iid, whole point
                                              hyper = interaction_hyper ))

ptm <- Sys.time()
typeI_fit <- inla(typeI_formula, data = ohio_df, family = 'poisson',
                  E = pop_at_risk, control.compute = list(config = TRUE, 
                                                          cpo = TRUE,    
                                                          waic = TRUE))
time_typeI = Sys.time() - ptm

#Inference on model w. type I interaction
print_cpo_etc(typeI_fit, time_typeI)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeI_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeI_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeI_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeI_fit)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeI_fit)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeI_fit)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeI_fit)

#Plot to check up on shit with the interactions to do ehh

#Precision of interactions
plot_prec_interactions(typeI_fit, 
                       "Posterior density of precision of type I interaction")

#Interactions themselves
plot_interaction(typeI_fit)


###

#Type II

#Get sum-to-zero constraints for type II interaction
typeII_constraints = constraints_maker(type = "II", n = n, t = T)

#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))
#Get precision matric for type II interaction by Kronecker product
typeII_prec <- scaled_RW_prec %x% diag(n)

typeII_formula <- update(base_formula, ~. + f(space_time_unstructured, 
                                              model = "generic0", 
                                              Cmatrix = typeII_prec, 
                                              extraconstr = typeII_constraints, 
                                              rankdef = n, 
                                              hyper = interaction_hyper))



ptm <- Sys.time()
typeII_fit <- inla(typeII_formula, data = ohio_df, family = "poisson",
                   E = pop_at_risk, control.compute = list(config = TRUE,
                                                           cpo = TRUE,
                                                           waic = TRUE))
time_typeII = Sys.time() - ptm

#Inference on type II model
print_cpo_etc(typeII_fit, time_typeII)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeII_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeII_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeII_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeII_fit)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeII_fit)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeII_fit)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeII_fit)


#Precision of interactions
plot_prec_interactions(typeII_fit, 
                       "Posterior density of precision of type II interaction")

#Interactions themselves
plot_interaction(typeII_fit)


###

#Type III

#Get constraints for the type III interactions
typeIII_constraints <- constraints_maker(type = "III", n = n, t = T)

# get scaled ICAR
scaled_ICAR_prec <- INLA::inla.scale.model(ICAR_prec, 
                     constr = list(A = matrix(1,1,dim(ICAR_prec)[1]), e = 0))

# Kronecker product between IID x ICAR
typeIII_prec <- diag(T) %x% scaled_ICAR_prec 

typeIII_formula <- update(base_formula, ~. + f(space_time_unstructured, 
                                               model = "generic0", 
                                               Cmatrix = typeIII_prec, 
                                               extraconstr = typeIII_constraints, 
                                               rankdef = T, 
                                               hyper = interaction_hyper))


ptm <- Sys.time()
typeIII_fit <- inla(typeIII_formula, data = ohio_df, family = "poisson",
                    E = pop_at_risk, control.compute = list(config = TRUE, 
                                                            cpo = TRUE,
                                                            waic = TRUE))
time_typeIII = Sys.time() - ptm

#Inference on type III model
print_cpo_etc(typeIII_fit, time_typeIII)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeIII_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeIII_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeIII_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeIII_fit)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeIII_fit)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeIII_fit)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeIII_fit)


#Precision of interactions
plot_prec_interactions(typeIII_fit, 
                       "Posterior density of precision of type III interaction")

#Interactions themselves
plot_interaction(typeIII_fit)

###

#Type IV

#Get constraints for type IV interactions
typeIV_constraints <- constraints_maker(type = "IV", n = n, t = T)

#Get type IV interaction precision matrix
typeIV_prec <- scaled_RW_prec %x% scaled_ICAR_prec

#Get formula for type IV
typeIV_formula <- update(base_formula, ~. + f(space_time_unstructured, 
                                              model = "generic0",
                                              Cmatrix = typeIV_prec,
                                              extraconstr = typeIV_constraints,
                                              rankdef = (n + T - 1), 
                                              hyper = interaction_hyper))


ptm <- Sys.time()
typeIV_fit <- inla(typeIV_formula, data = ohio_df, family = "poisson",
                    E = pop_at_risk, control.compute = list(config = TRUE, 
                                                            cpo = TRUE,
                                                            waic = TRUE))
time_typeIV = Sys.time() - ptm

#Inference on type IV model

print_cpo_etc(typeIV_fit, time_typeIV)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeIV_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeIV_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeIV_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeIV_fit)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeIV_fit)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeIV_fit)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeIV_fit)


#Precision of interactions
plot_prec_interactions(typeIV_fit, 
                       "Posterior density of precision of type IV interaction")

#Interactions themselves
plot_interaction(typeIV_fit)

