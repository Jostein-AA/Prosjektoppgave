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

#Extract adjacency structure and create adjacency matrix
nb <- spdep::poly2nb(ohio_map, queen = FALSE)
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
adjacency_matrix <- Matrix(matrix4inla, sparse = TRUE) #Make it sparse

#get RW1 prec
#Get precision matricies for RW1
RW1_prec <- INLA:::inla.rw(n = T, order = 1, 
                           scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                           sparse = TRUE)

#Take Kronecker product to get graph structure

#Scale precision matrix of RW model so the geometric mean of the marginal variances is one
scaled_RW_prec <- inla.scale.model(RW1_prec,
                                   list(A = matrix(1, 1, dim(RW1_prec)[1]),
                                        e = 0))

# get scaled ICAR
scaled_adjacency_matrix <- INLA::inla.scale.model(adjacency_matrix, 
                                                  constr = list(A = matrix(1,1,dim(adjacency_matrix)[1]), e = 0))

x_graph <- scaled_RW_prec %x% scaled_adjacency_matrix 

epsilon_graph <- diag(T) %x% scaled_adjacency_matrix 

#Define prior
alt_hyper = list(prec = list(prior="pc.prec",
                             param=c(1,0.01)),
                 lambda = list(prior = "gaussian",
                               param = c(0, 0.5)))

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

#This formulation is third best model, not needing any constraints (except random effects).
#Analouge type IV interaction though! Not really what I want in the same way as what has been done
alternate_formula <- deaths ~ 1 + f(year, 
                                    model = 'bym2',
                                    scale.model = T, 
                                    constr = T, 
                                    rankdef = 1,
                                    graph = scaled_RW_prec,
                                    hyper = temporal_hyper) + 
                                  f(county, 
                                    model = 'bym2',
                                    scale.model = T,
                                    constr = T,
                                    rankdef = 1,
                                    graph = scaled_adjacency_matrix,
                                    hyper = spatial_hyper) + 
                                  f(space.time, 
                                    model = "besagproper2",
                                    graph = x_graph,
                                    hyper = alt_hyper) 


#Fit the base formula 
ptm <- Sys.time()
alternate_fit <- inla(alternate_formula,
                      data = ohio_df,
                      family = "poisson",
                     E = pop_at_risk, 
                     control.compute = list(config = TRUE, # To see constraints later
                                            cpo = T,   # For model selection
                                            waic = T)) # For model selection

time_alt = Sys.time()-ptm
print(c("Alternate model fitted in: ", time_alt))
print(mean(-log(alternate_fit$cpo$cpo)))




