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
  RW2_prec <- INLA:::inla.rw(n = time, order = 2, 
                             scale.model = FALSE, 
                             sparse = TRUE)
  
  #Make the base formula
  base_formula <- deaths ~ 1 + f(year, 
                                 model = 'bym2',
                                 scale.model = T, 
                                 constr = T, 
                                 rankdef = 2,
                                 graph = RW2_prec,
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
  
  
  RW2_ICAR_I_fit <- inla(typeI_formula,
                         data = temp_ohio,
                         family = 'poisson',
                         E = pop_at_risk,
                         control.predictor = list(compute = TRUE),       #For predictions
                         control.compute = list(config = TRUE, # To see constraints later
                                                cpo = T,   # For model selection
                                                waic = T,  # For model selection
                                                return.marginals.predictor=TRUE)) #For predictions
  print("I")
  
  
  
  #Scale precision matrix of RW model so the geometric mean of the marginal variances is one
  scaled_RW_prec <- inla.scale.model(RW2_prec,
                                     list(A = matrix(1, 1, dim(RW2_prec)[1]),
                                          e = 0))
  #Get precision matric for type II interaction by Kronecker product
  typeII_prec <- scaled_RW_prec %x% diag(n)
  
  
  #Get constraints for type II interaction using a RW2
  typeII_constraints = constraints_maker(type = "II", n = n, t = time, 
                                         rw = "RW2", prec_matrix = typeII_prec)
  
  
  
  typeII_formula <- update(base_formula, ~. + f(space.time, 
                                                model = "generic0", 
                                                Cmatrix = typeII_prec, 
                                                extraconstr = typeII_constraints, 
                                                rankdef = (2*n), 
                                                hyper = interaction_hyper))
  
  
  RW2_ICAR_II_fit <- inla(typeII_formula,
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
  
  
  
  RW2_ICAR_III_fit <- inla(typeIII_formula,
                           data = temp_ohio,
                           family = "poisson",
                           control.predictor = list(compute = TRUE),       #For predictions
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T,   # For model selection
                                                  waic = T,  # For model selection
                                                  return.marginals.predictor=TRUE)) #For predictions
  
  print("III")
  
  #Get type IV interaction precision matrix
  typeIV_prec <- scaled_RW_prec %x% scaled_ICAR_prec
  
  
  typeIV_constraints <- constraints_maker(type = "IV", n = n, t = time,
                                          rw = "RW2", prec_matrix = typeIV_prec)
  
  
  
  #Get formula for type IV
  typeIV_formula <- update(base_formula, ~. + f(space.time, 
                                                model = "generic0",
                                                Cmatrix = typeIV_prec,
                                                extraconstr = typeIV_constraints,
                                                rankdef = (2 * n + time - 2), 
                                                hyper = interaction_hyper))
  
  
  RW2_ICAR_IV_fit <- inla(typeIV_formula,
                          data = temp_ohio,
                          family = "poisson",
                          E = pop_at_risk,
                          control.predictor = list(compute = TRUE),       #For predictions
                          control.compute = list(config = TRUE, # To see constraints later
                                                 cpo = T,   # For model selection
                                                 waic = T,  # For model selection
                                                 return.marginals.predictor=TRUE)) #For predictions
  
  print("IV")
  
  
  #Extract marginal of linear predictor predictive distribution
  #Save the linear predictor marginals because this is what we will use for CRPS
  base_prediction_marginals = lapply(base_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal) #base_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)] 
  I_prediction_marginals =    lapply(RW2_ICAR_I_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_I_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  II_prediction_marginals =   lapply(RW2_ICAR_II_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_II_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  III_prediction_marginals =  lapply(RW2_ICAR_III_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_III_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  IV_prediction_marginals =   lapply(RW2_ICAR_IV_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                     FUN = my_inla_t_marginal)#RW1_ICAR_IV_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)]
  
  if(time == 12){
    improper_2_noInt_pred = base_prediction_marginals
    improper_2_I_pred    = I_prediction_marginals
    improper_2_II_pred   = II_prediction_marginals
    improper_2_III_pred  = III_prediction_marginals
    improper_2_IV_pred   = IV_prediction_marginals
    
    
  } else {
    improper_2_noInt_pred = rbind(improper_2_noInt_pred,
                                  base_prediction_marginals)
    
    improper_2_I_pred = rbind(improper_2_I_pred,
                        I_prediction_marginals)
    
    improper_2_II_pred = rbind(improper_2_II_pred,
                         II_prediction_marginals)
    
    improper_2_III_pred = rbind(improper_2_III_pred,
                          III_prediction_marginals)
    
    improper_2_IV_pred = rbind(improper_2_IV_pred,
                         IV_prediction_marginals)
  }
  
  print(paste("Predicted for time: ", time, "/ Time used so far: ", Sys.time() - ptm))
}
print(paste("Time used on predicting improper models: ", Sys.time() - ptm))


save(improper_2_noInt_pred,
     improper_2_I_pred,
     improper_2_II_pred,
     improper_2_III_pred,
     improper_2_IV_pred,
     file = "one_step_predictions_RW2.RData")

