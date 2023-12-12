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


#Function used because specifying link function did not work for me, 
#for some reason. 
my_inla_t_marginal <- function(prediction_marginal){
  return(inla.tmarginal(function(x){exp(x)}, prediction_marginal))
}


for(time in 12:21){ #For loop to sequentially predict one and one year ahead, start at year = 12, and go until end
  temp_ohio = ohio_df[ohio_df$year<= time, ] #Extract data in year 1:time
  temp_ohio[temp_ohio$year == time, ]$deaths = NA
  
  #Get constraints for the type III interactions
  A <- matrix(0, nrow = time, ncol = n * time)
  for (i in 1:time) {
    # The ICAR at each time point needs to sum to 0
    A[i, ((i - 1) * n + 1):(i * n)] <- 1
  }
  
  typeIII_constraints <- list(A = A, e = rep(0, dim(A)[1]))
  
  
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
                                                  return.marginals.predictor=TRUE), #For predictions
                           verbose = T)
  
  print("1 III")
  
  RW2_prec <- INLA:::inla.rw(n = time, order = 2, 
                             scale.model = FALSE, 
                             sparse = TRUE)
  
  
  #Make the base formula
  base_formula_2 <- deaths ~ 1 + f(year, 
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
  
  
  RW2_typeIII_formula <- update(base_formula_2, ~. + f(space.time, 
                                                       model = "generic0", 
                                                       Cmatrix = typeIII_prec, 
                                                       extraconstr = typeIII_constraints, 
                                                       rankdef = time, 
                                                       hyper = interaction_hyper))
  
  RW2_ICAR_III_fit <- inla(RW2_typeIII_formula,
                           data = temp_ohio,
                           family = "poisson",
                           control.predictor = list(compute = TRUE),       #For predictions
                           control.compute = list(config = TRUE, # To see constraints later
                                                  cpo = T,   # For model selection
                                                  waic = T,  # For model selection
                                                  return.marginals.predictor=TRUE), #For predictions
                           verbose = T)

  print("2 III")
  
  #Extract marginal of linear predictor at times predicted on,
  #use inla.tmarginal to transform to exp(linear predictor)
  RW1_III_prediction_marginals =  lapply(RW1_ICAR_III_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                   FUN = my_inla_t_marginal)
  
  RW2_III_prediction_marginals =  lapply(RW2_ICAR_III_fit$marginals.fitted.values[(n * (time - 1) + 1):(n * time)],
                                         FUN = my_inla_t_marginal)
 
  if(time == 12){
    RW1_III_predicted  = RW1_III_prediction_marginals
    RW2_III_predicted  = RW2_III_prediction_marginals
    
    
  } else {
    RW1_III_predicted = rbind(RW1_III_predicted,
                          RW1_III_prediction_marginals)
    
    RW2_III_predicted = rbind(RW2_III_predicted,
                              RW2_III_prediction_marginals)
    
  }
  
  print(paste("Predicted for time: ", time))
}












