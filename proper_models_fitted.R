#Clear environment
rm(list = ls())

#Load necessary libraries
library(INLA)
library(spData)
library(sf)
library(spdep)

#Load utility functions
#source("utilities.R")

#Load data
ohio_df <- read.csv("ohio_df.csv")
ohio_df$year <- ohio_df$year - min(ohio_df$year) + 1

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]

#Get number of counties (n) and number of years (T)
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points


#Make precision matrix for AR(1)
a = 0.99 #Close to 1, but will still result in proper distributions
ar1_prec = diag(1+a**2, nrow = T)
ar1_prec[1, 1] = 1; ar1_prec[T, T] = 1
ar1_prec[1, 2] = -a; ar1_prec[T, T - 1] = -a
for(i in 2:(T-1)){
  ar1_prec[i, i + 1] = -a; ar1_prec[i, i - 1] = -a
}
ar1_prec = Matrix(ar1_prec, sparse = TRUE) #Make it sparse

#Make precision matrix for CAR (We will have same as ICAR but ensure that it is Diagonal Dominant)
nb <- spdep::poly2nb(ohio_map, queen = FALSE)
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla) + 0.5
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
CAR_prec <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse


#Specify hyperparameters with corresponding priors
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



#Make base formula
#Make the base formula
#base_formula <- deaths ~ 1 + year +
#                             f(year, 
#                               model = 'bym2',
#                               scale.model = T, 
#                               constr = F, 
#                               rankdef = 0,
#                               graph = ar1_prec,
#                               hyper = temporal_hyper) + 
#                             f(county, 
#                               model = 'bym2',
#                               scale.model = T,
#                               constr = F,
#                               rankdef = 0,
#                               graph = CAR_prec,
#                               hyper = spatial_hyper)


base_formula <- deaths ~ 1 + f(year, model = "ar1") +
                             f(county, model = "besagproper2")

###

#Fit the base formula 
ptm <- Sys.time()
ar1_CAR_fit <- inla(base_formula, data = ohio_df, family = "poisson",
                     E = pop_at_risk, 
                     control.compute = list(config = TRUE, # To see constraints later
                                            cpo = T,   # For model selection
                                            waic = T)) # For model selection

time_AR1_CAR = Sys.time()-ptm
print(c("Basic model fitted in: ", time_AR1_CAR))

###

#Update base formula to also contain iid interaction
typeI_formula <- update(base_formula,  ~. + f(space.time,
                                              model="iid", #Has to be iid, whole point
                                              hyper = interaction_hyper ))

ptm <- Sys.time()
RW1_ICAR_I_fit <- inla(typeI_formula, data = ohio_df, family = 'poisson',
                       E = pop_at_risk, control.compute = list(config = TRUE, 
                                                               cpo = TRUE,    
                                                               waic = TRUE))
time_RW1_ICAR_I = Sys.time() - ptm
print(c("Type I model fitted in: ", time_RW1_ICAR_I))


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

typeII_formula <- update(base_formula, ~. + f(space.time, 
                                              model = "generic0", 
                                              Cmatrix = typeII_prec, 
                                              extraconstr = typeII_constraints, 
                                              rankdef = n, 
                                              hyper = interaction_hyper))



ptm <- Sys.time()
RW1_ICAR_II_fit <- inla(typeII_formula, data = ohio_df, family = "poisson",
                        E = pop_at_risk, control.compute = list(config = TRUE,
                                                                cpo = TRUE,
                                                                waic = TRUE))
time_RW1_ICAR_II = Sys.time() - ptm
print(c("Type II model fitted in: ", time_RW1_ICAR_II))


###

#Type III

#Get constraints for the type III interactions
typeIII_constraints <- constraints_maker(type = "III", n = n, t = T)

# get scaled ICAR
scaled_ICAR_prec <- INLA::inla.scale.model(ICAR_prec, 
                                           constr = list(A = matrix(1,1,dim(ICAR_prec)[1]), e = 0))

# Kronecker product between IID x ICAR
typeIII_prec <- diag(T) %x% scaled_ICAR_prec 

typeIII_formula <- update(base_formula, ~. + f(space.time, 
                                               model = "generic0", 
                                               Cmatrix = typeIII_prec, 
                                               extraconstr = typeIII_constraints, 
                                               rankdef = T, 
                                               hyper = interaction_hyper))



ptm <- Sys.time()
RW1_ICAR_III_fit <- inla(typeIII_formula, data = ohio_df, family = "poisson",
                         E = pop_at_risk, control.compute = list(config = TRUE, 
                                                                 cpo = TRUE,
                                                                 waic = TRUE))

time_RW1_ICAR_III = Sys.time() - ptm
print(c("Type III model fitted in: ", time_RW1_ICAR_III))


###

#Type IV

#Get constraints for type IV interactions
typeIV_constraints <- constraints_maker(type = "IV", n = n, t = T)

#Get type IV interaction precision matrix
typeIV_prec <- scaled_RW_prec %x% scaled_ICAR_prec

#Get formula for type IV
typeIV_formula <- update(base_formula, ~. + f(space.time, 
                                              model = "generic0",
                                              Cmatrix = typeIV_prec,
                                              extraconstr = typeIV_constraints,
                                              rankdef = (n + T - 1), 
                                              hyper = interaction_hyper))


ptm <- Sys.time()
RW1_ICAR_IV_fit <- inla(typeIV_formula, data = ohio_df, family = "poisson",
                        E = pop_at_risk, control.compute = list(config = TRUE, 
                                                                cpo = TRUE,
                                                                waic = TRUE))
time_RW1_ICAR_IV = Sys.time() - ptm
print(c("Type IV model fitted in: ", time_RW1_ICAR_IV))



#Save INLA objects
save(n, T, ohio_map, ohio_df,
     RW1_ICAR_fit, time_RW1_ICAR,
     RW1_ICAR_I_fit, time_RW1_ICAR_I,
     RW1_ICAR_II_fit, time_RW1_ICAR_II,
     RW1_ICAR_III_fit, time_RW1_ICAR_III,
     RW1_ICAR_IV_fit, time_RW1_ICAR_IV,
     file = "improper_RW1_ICAR_fitted.RData")





