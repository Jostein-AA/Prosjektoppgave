#Clear environment
rm(list = ls())

#Load necessary libraries
library(INLA)
library(spData)
library(sf)
library(spdep)

#Load utility functions
source("utilities.R")

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

####

#Extract adjacency structure and create precision matrix for ICAR
nb <- spdep::poly2nb(ohio_map, queen = FALSE)
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
ICAR_prec <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse


# get scaled ICAR
scaled_ICAR_prec <- INLA::inla.scale.model(ICAR_prec, 
                                           constr = list(A = matrix(1,1,dim(ICAR_prec)[1]), e = 0))


#Define the different formulas 
#Temporal hyperparameters (Precision of iid and precision of RW1) w. corresponding priors: penalized constraint 
temporal_hyper = list(prec = list(prior = 'pc.prec', 
                                  param = c(1, 0.01)), #Magic numbers
                      phi = list(prior = 'pc', 
                                 param = c(0.5, 0.5)) #Magic numbers
) 

#Spatial hyperparameters (Precision of iid and precision of ICAR) w. corresponding priors: penalized constraint
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.01)), #Magic numbers
                     phi = list(prior = 'pc', 
                                param = c(0.5, 0.5)) #Magic numbers
)
#Interaction hyperparameter and prior (Precision of interaction)
interaction_hyper = list(theta=list(prior="pc.prec",
                                    param=c(1,0.01)))



####
#Improper models fitted
ptm <- Sys.time()
for(time in 12:21){ #For loop to sequentially predict one and one year ahead, start at year = 11, and go until end
  temp_ohio = ohio_df[ohio_df$year<= time, ] #Extract data in year 1:time
  temp_ohio[temp_ohio$year == time, ]$deaths = NA
  
  #make RW1 prec for these years
  RW1_prec <- INLA:::inla.rw(n = time, order = 1, 
                             scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                             sparse = TRUE)
  
  #Make the base formula
  base_formula <- deaths ~ 1 + f(year, 
                                 model = 'bym2',
                                 scale.model = T, 
                                 constr = T, 
                                 rankdef = 1,
                                 graph = RW1_prec,
                                 hyper = temporal_hyper) + 
                               f(county, 
                                 model = 'bym2',
                                 scale.model = T,
                                 constr = T,
                                 rankdef = 1,
                                 graph = ICAR_prec,
                                 hyper = spatial_hyper)
  
  base_fit <- inla(base_formula, 
                   data = temp_ohio,
                   family = "poisson",
                   E = pop_at_risk, 
                   control.predictor = list(compute = TRUE),       #For predictions
                   control.compute = list(config = TRUE, # To see constraints later
                                          cpo = T,   # For model selection
                                          waic = T,  # For model selection
                                          return.marginals.predictor=TRUE)) #For predictions
  
  print("base")
  
  
  #Update base formula to also contain iid interaction
  typeI_formula <- update(base_formula,  ~. + f(space.time,
                                                model="iid", #Has to be iid, whole point
                                                hyper = interaction_hyper ))
  
  
  RW1_ICAR_I_fit <- inla(typeI_formula,
                         data = temp_ohio,
                         family = 'poisson',
                         E = pop_at_risk,
                         control.predictor = list(compute = TRUE),       #For predictions
                         control.compute = list(config = TRUE, # To see constraints later
                                                cpo = T,   # For model selection
                                                waic = T,  # For model selection
                                                return.marginals.predictor=TRUE)) #For predictions
  print("I")
  
  typeII_constraints = constraints_maker(type = "II", n = n, t = time)
  
  #Scale precision matrix of RW model so the geometric mean of the marginal variances is one
  scaled_RW_prec <- inla.scale.model(RW1_prec,
                                     list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                          e = 0))
  #Get precision matric for type II interaction by Kronecker product
  typeII_prec <- scaled_RW_prec %x% diag(n)
  
  typeII_formula <- update(base_formula, ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeII_prec, 
                                                extraconstr = typeII_constraints, 
                                                rankdef = n, 
                                                hyper = interaction_hyper))
  
  
  RW1_ICAR_II_fit <- inla(typeII_formula,
                          data = temp_ohio,
                          family = "poisson",
                          E = pop_at_risk, 
                          control.predictor = list(compute = TRUE),       #For predictions
                          control.compute = list(config = TRUE, # To see constraints later
                                                 cpo = T,   # For model selection
                                                 waic = T,  # For model selection
                                                 return.marginals.predictor=TRUE)) #For predictions
  print("II")
  
  
  #Get constraints for the type III interactions
  typeIII_constraints <- constraints_maker(type = "III", n = n, t = time)
  
  
  # Kronecker product between IID x ICAR
  typeIII_prec <- diag(time) %x% scaled_ICAR_prec 
  
  typeIII_formula <- update(base_formula, ~. + f(space.time, 
                                                 model = "generic0", 
                                                 Cmatrix = typeIII_prec, 
                                                 extraconstr = typeIII_constraints, 
                                                 rankdef = time, 
                                                 hyper = interaction_hyper))
  
  
  
  RW1_ICAR_III_fit <- inla(typeIII_formula,
                           data = temp_ohio,
                           family = "poisson",
                           control.predictor = list(compute = TRUE),       #For predictions
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T,   # For model selection
                                                  waic = T,  # For model selection
                                                  return.marginals.predictor=TRUE)) #For predictions
  
  print("III")
  
  typeIV_constraints <- constraints_maker(type = "IV", n = n, t = time)
  
  #Get type IV interaction precision matrix
  typeIV_prec <- scaled_RW_prec %x% scaled_ICAR_prec
  
  #Get formula for type IV
  typeIV_formula <- update(base_formula, ~. + f(space.time, 
                                                model = "generic0",
                                                Cmatrix = typeIV_prec,
                                                extraconstr = typeIV_constraints,
                                                rankdef = (n + time - 1), 
                                                hyper = interaction_hyper))
  
  
  RW1_ICAR_IV_fit <- inla(typeIV_formula,
                          data = temp_ohio,
                          family = "poisson",
                          E = pop_at_risk,
                          control.predictor = list(compute = TRUE),       #For predictions
                          control.compute = list(config = TRUE, # To see constraints later
                                                 cpo = T,   # For model selection
                                                 waic = T,  # For model selection
                                                 return.marginals.predictor=TRUE)) #For predictions
  
  print("IV")
  
  
  #base_predicted = lapply(base_prediction_marginals, FUN = my_inla_t_marginal)
  #I_predicted = lapply(I_prediction_marginals, FUN = my_inla_t_marginal)
  #II_predicted = lapply(II_prediction_marginals, FUN = my_inla_t_marginal)
  #III_predicted = lapply(III_prediction_marginals, FUN = my_inla_t_marginal)
  #IV_predicted = lapply(IV_prediction_marginals, FUN = my_inla_t_marginal)
  
  #Extract marginal of linear predictor predictive distribution
  #Save the linear predictor marginals because this is what we will use for CRPS
  base_prediction_marginals = lapply(base_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal) #base_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)] 
  I_prediction_marginals =    lapply(RW1_ICAR_I_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_I_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  II_prediction_marginals =   lapply(RW1_ICAR_II_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_II_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  III_prediction_marginals =  lapply(RW1_ICAR_III_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_III_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  IV_prediction_marginals =   lapply(RW1_ICAR_IV_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_IV_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  
  if(time == 12){
    base_predicted = base_prediction_marginals
    I_predicted    = I_prediction_marginals
    II_predicted   = II_prediction_marginals
    III_predicted  = III_prediction_marginals
    IV_predicted   = IV_prediction_marginals
    
    
  } else {
    base_predicted = rbind(base_predicted,
                           base_prediction_marginals)
    
    I_predicted = rbind(I_predicted,
                        I_prediction_marginals)
    
    II_predicted = rbind(II_predicted,
                         II_prediction_marginals)
    
    III_predicted = rbind(III_predicted,
                          III_prediction_marginals)
    
    IV_predicted = rbind(IV_predicted,
                         IV_prediction_marginals)
  }
  
  print(paste("Predicted for time: ", time, "/ Time used so far: ", Sys.time() - ptm))
}
print(paste("Time used on predicting improper models: ", Sys.time() - ptm))


####
#Proper models now
ohio_df_changed <- ohio_df

#Firstly, switch ordering from year over county, to county over year
#This is done for the sake of the space time interaction in this model
ohio_df_changed = ohio_df_changed[order(ohio_df_changed$county, decreasing = F), ] 
ohio_df_changed$space.time <- 1:(n *T); rownames(ohio_df_changed) <- 1:nrow(ohio_df_changed) 

#Make copies of county and year.
#It is done because both year and county is used twice in some formulas
ohio_df_changed$county.copy <- ohio_df_changed$county
ohio_df_changed$year.copy <- ohio_df_changed$year


#Define hyperparameters and corresponding priors
#Define Temporal hyperparameters and corresponding priors
ar1_hyper = list(prec = list(prior = 'pc.prec', 
                             param = c(1, 0.008)), #Magic numbers
                 rho = list(prior = 'normal', 
                            param = c(0, 0.5)), #Magic numbers
                 mean = list(prior = 'normal',
                             param = c(0, 2.5))) #Magic numbers


#Define Spatial hyperparameters and corresponding priors
spatial_hyper = list(prec= list(prior = 'pc.prec', 
                                param = c(1, 0.008)), #Magic numbers
                     lambda = list(prior = 'gaussian', 
                                   param = c(0, 1))) #Magic numbers


#Define Temporal hyperparameters and corresponding priors
ar1_hyper_full = list(prec = list(prior = 'pc.prec', 
                                  param = c(1, 0.005)), #Magic numbers
                      rho = list(prior = 'normal', 
                                 param = c(0, 1)), #Magic numbers, second is precision
                      mean = list(prior = 'normal',
                                  param = c(0, 2))) #Magic numbers


#Define Spatial hyperparameters and corresponding priors
spatial_hyper_full = list(prec= list(prior = 'pc.prec', 
                                    param = c(1, 0.005)), #Magic numbers
                          lambda = list(prior = 'gaussian', 
                                        param = c(0, 1))) #Magic numbers



ptm <- Sys.time()
for(time in 12:21){ #For loop to sequentially predict one and one year ahead, start at year = 11, and go until end
  temp_ohio = ohio_df_changed[ohio_df_changed$year <= time, ] #Extract data in year 1:time
  temp_ohio[temp_ohio$year == time, ]$deaths = NA; rownames(temp_ohio) <- 1:nrow(temp_ohio)
  
  
  proper_base_formula <- deaths ~ 1 + year +
                                  f(year.copy,
                                    model = "ar1",
                                    hyper = ar1_hyper) + 
                                  f(county, 
                                    model = "besagproper2",
                                    graph = ICAR_prec,
                                    hyper = spatial_hyper)
  
  proper_base <- inla(proper_base_formula,
                      data = temp_ohio,
                      family = "poisson",
                      E = pop_at_risk,
                      control.predictor = list(compute = TRUE),       #For predictions
                      control.compute = list(config = TRUE, # To see constraints later
                                             cpo = T,   # For model selection
                                             waic = T,  # For model selection
                                             return.marginals.predictor=TRUE)) #For predictions
  
  print("proper base")
  
  proper_interaction_formula <- deaths ~ 1 + year + 
                                         f(county, 
                                           model = "besagproper2",
                                           graph = ICAR_prec,
                                           hyper = spatial_hyper,
                                           group = year, 
                                           control.group = list(model = "ar1"))
  
  proper_interaction <- inla(proper_interaction_formula,
                             data = temp_ohio,
                             family = "poisson",
                             E = pop_at_risk,
                             control.predictor = list(compute = TRUE),       #For predictions
                             control.compute = list(config = TRUE, # To see constraints later
                                                    cpo = T,   # For model selection
                                                    waic = T,  # For model selection
                                                    return.marginals.predictor=TRUE)) #For predictions
  
  print("proper interaction")
  
  proper_full_formula <- deaths ~ 1 + year + 
                                  f(year.copy,
                                    model = "ar1",
                                    hyper = ar1_hyper_full) +
                                  f(county, 
                                    model = "besagproper2",
                                    graph = ICAR_prec,
                                    hyper = spatial_hyper_full) + 
                                  f(county.copy, 
                                    model = "besagproper2",
                                    graph = ICAR_prec,
                                    hyper = spatial_hyper_full,
                                    group = year, 
                                    control.group = list(model = "ar1")) 
  
  
  proper_full <- inla(proper_full_formula,
                      data = temp_ohio,
                      family = "poisson",
                      E = pop_at_risk,
                      control.predictor = list(compute = TRUE),       #For predictions
                      control.inla = list(int.strategy = "grid"),
                      control.compute = list(config = TRUE, # To see constraints later
                                             cpo = T,   # For model selection
                                             waic = T,  # For model selection
                                             return.marginals.predictor=TRUE)) #For predictions
  
  print("proper full")
  
  #Have to sort the fitted/predicted values because they are in different order to the improper ones
  proper_base_marginals         = sort_proper_fitted_2(proper_base$marginals.fitted.values, n, time)
  proper_interactions_marginals = sort_proper_fitted_2(proper_interaction$marginals.fitted.values, n, time)
  proper_full_marginals         = sort_proper_fitted_2(proper_full$marginals.fitted.values, n, time)
  
  #Save the linear predictor marginals because this is what we will use for CRPS, AE
  proper_base_marginals         = lapply(proper_base_marginals[(n * (time - 1) + 1):(n * time)],
                                         FUN = my_inla_t_marginal)#proper_base_marginals[(n * (time - 1) + 1):(n * time)]
  proper_interactions_marginals = lapply(proper_interactions_marginals[(n * (time - 1) + 1):(n * time)],
                                         FUN = my_inla_t_marginal)#proper_interactions_marginals[(n * (time - 1) + 1):(n * time)]
  proper_full_marginals         = lapply(proper_full_marginals[(n * (time - 1) + 1):(n * time)],
                                         FUN = my_inla_t_marginal)#proper_full_marginals[(n * (time - 1) + 1):(n * time)]
  
  if(time == 12){
    proper_base_predicted        = proper_base_marginals
    proper_interaction_predicted = proper_interactions_marginals
    proper_full_predicted        = proper_full_marginals
    
  } else {
    proper_base_predicted = rbind(proper_base_predicted,
                                  proper_base_marginals)
    proper_interaction_predicted = rbind(proper_interaction_predicted,
                                         proper_interactions_marginals)
    proper_full_predicted = rbind(proper_full_predicted,
                                  proper_full_marginals)
  }
  print(paste("Predicted for time: ", time, "/ Time used so far: ", Sys.time() - ptm))
}
print(paste("Time used on predicting proper models: ", Sys.time() - ptm))


####
#Save results

save(n, T, ohio_map, ohio_df, ohio_df_changed,
     base_predicted,
     I_predicted,
     II_predicted,
     III_predicted,
     IV_predicted,
     proper_base_predicted,
     proper_interaction_predicted,
     proper_full_predicted,
     file = "one_step_predictions.RData")



























