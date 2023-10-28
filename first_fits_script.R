#Load libraries
library(INLA)
library(tidyverse)
library(spData)
library(sf)
library(spdep)
library(ggplot2)
library(ggspatial)
require(mgcv)
library(MASS)
library(fastmatrix)
library(rwc)
library(gridExtra)
library(ggpubr)
library(roahd)
###

#Load data
###
#Load structure matrices
struct_RW1 = readMM(file='struct_RW1.txt')
struct_RW2 = readMM(file = 'struct_RW2.txt')
ICAR_structure = readMM(file = 'ICAR_struct.txt')

#Load data
ohio_df <- read.csv("ohio_df.csv")

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]

#Calculate adjacencies (if needed)
nb <- spdep::poly2nb(ohio_map, queen = FALSE)

#Get number of counties (n) and number of years (T)
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points
###

#Specify hyperparameter priors
###
#Temporal
temporal_hyper = list(prec.unstruct = list(prior = 'pc.prec', 
                                           param = c(1, 0.01)), #Magic numbers
                      prec.spatial = list(prior = 'pc.prec', 
                                          param = c(1, 0.01)) #Magic numbers
) 

#Spatial
spatial_hyper = list(prec.unstruct = list(prior = 'pc.prec', 
                                          param = c(1, 0.01)), #Magic numbers
                     prec.spatial = list(prior = 'pc.prec', 
                                         param = c(1, 0.01)) #Magic numbers
)
###

#Basic formula
###
#Make the base formula (is the overall mean specified correctly here???)
basic_linear_formula <- deaths ~ 1 + f(year, 
                                       model = 'bym',
                                       scale.model = T, 
                                       constr = T, 
                                       rankdef = 1,
                                       graph = struct_RW1,
                                       hyper = temporal_hyper
                                       ) + 
                                     f(county, 
                                       model = 'bym',
                                       scale.model = T,
                                       constr = T,
                                       rankdef = 1,
                                       graph = ICAR_structure,
                                       hyper = spatial_hyper
                                       )
###

#Fit the model w.o. interactions
###
#Measure time to do inference
ptm <- Sys.time()
basic_model_fit <- inla(basic_linear_formula,
                        data = ohio_df,
                        family = "poisson",
                        E = pop_at_risk, 
                        control.compute = list(config = TRUE, # needed if you want to see constraints later
                                               cpo = T,
                                               waic = T), #Needed for model selection
                        verbose = F #Not T unless want to see a lot of stuff
) 
time = Sys.time()-ptm

#Check that fitted values make sense, check random effects. Check that time effect makes sense.

print(c("Number of constraints (should be 2): ",toString(basic_model_fit$misc$configs$constr$nc)))
print(c("mean CPO:", toString(round(mean(basic_model_fit$cpo$cpo), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(basic_model_fit$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(basic_model_fit$cpo$failure[1])))
print(c("CPU: ", toString(summary(basic_model_fit)$cpu.used)))
print(c("Time to compute", toString(time)))



#Plot just the temporal effects
plot(basic_model_fit$summary.random$year$ID[1:T], 
     basic_model_fit$summary.random$year$mean[1:T], 
     type = "l", lty = 1, xlab = "Year: t", ylab = "random effect alpha_t") + 
  lines(basic_model_fit$summary.random$year$ID[1:T],
        basic_model_fit$summary.random$year$'0.025quant'[1:T],
        type = "l", lty = 2, col = "red") + 
  lines(basic_model_fit$summary.random$year$ID[1:T],
        basic_model_fit$summary.random$year$'0.975quant'[1:T],
        type = "l", lty = 2, col = "green") 

#Plot the spatial effects
spatial_structured_effect_mean <- basic_model_fit$summary.random$county$mean[1:n]
spatial_structured_effect_q025 <- basic_model_fit$summary.random$county$'0.025quant'[1:n]
spatial_structured_effect_q975 <- basic_model_fit$summary.random$county$'0.975quant'[1:n]
spatial_structured_effect_sd   <- basic_model_fit$summary.random$county$sd[1:n]

temp_ohio_map <- ohio_map
temp_ohio_map$mean <- spatial_structured_effect_mean
temp_ohio_map$q025 <- spatial_structured_effect_q025
temp_ohio_map$q975 <- spatial_structured_effect_q975
temp_ohio_map$sd   <- spatial_structured_effect_sd

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red

ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = mean), 
          alpha = 1,
          color="black") + ggtitle("spatial structured effect mean each county") +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_blank(), #Remove axis and background grid
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin =  unit(c(0, 0, 0, 0), "inches"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.spacing = unit(1, 'lines')) +
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
  binned_scale( #Scaling the color
    aesthetics = "fill",
    scale_name = "gradientn",
    palette = function(x) c(scale),
    labels = function(x){x},
    guide = "colorscale")

ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = sd), 
          alpha = 1,
          color="black") + ggtitle("spatial structured effect standard devation each county") +
  theme(plot.title = element_text(size = 20),
        axis.title.x = element_blank(), #Remove axis and background grid
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin =  unit(c(0, 0, 0, 0), "inches"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
        panel.spacing = unit(1, 'lines')) +
  guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
  binned_scale( #Scaling the color
    aesthetics = "fill",
    scale_name = "gradientn",
    palette = function(x) c(scale),
    labels = function(x){x},
    guide = "colorscale")

#Look at the fitted values compared to the actual values
print(basic_model_fit$summary.fitted.values)


###
#Fit model w. type I interactions
# Type I Interaction:
interaction_hypers = list(prec = list(param = c(1, 0.01)))

typeI_formula <- update(basic_linear_formula, 
                        ~. + f(space_time_unstructured,
                               model="iid", #Has to be iid, whole point
                               hyper = interaction_hypers ))





ptm <- Sys.time()
typeI_fit <- inla(typeI_formula,
                  data = ohio_df,
                  family = 'poisson',
                  E = pop_at_risk,
                  control.compute = list(config = TRUE, # needed if you want to see constraints later
                                         cpo = TRUE,    #Needed for model choice
                                         waic = TRUE), 
                  verbose = F #Not T unless want to see a lot of stuff
)
time = Sys.time() - ptm
print(c("Number of constraints (should be 2): ",toString(typeI_fit$misc$configs$constr$nc)))
print(c("mean CPO:", toString(round(mean(typeI_fit$cpo$cpo), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(typeI_fit$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(typeI_fit$cpo$failure[1])))
print(c("CPU: ", toString(summary(typeI_fit)$cpu.used)))
print(c("Time: ", time))
###

#Fit model w. type II interaction
###
# Type II: RW-time x iid-space
# creating constraint matrix needed for Type II Interaction

#sum-to-zero constraint be wilding (Is it the right way??? Does it actually do what it is supposed to???)
# - n rows as this interaction makes it so there is a RW for each county independent of each other
# - n * T columns as there are a total of n * T interaction terms (delta_(it))
# - sets A[i, which((0:(n * T - 1))%%n == i - 1)] = 1, as for county i, the interaction (delta_(it))
#     is only dependent on (delta_(i,t+-1)). Therefore the sum-to-zero is over these 88 RWs
# - Constraints: The RW1 in each area needs to sum to 0. Hence in constr.st e is a zero vector
A <- matrix(0, nrow = n, ncol = n * T)
for (i in 1:n) {
  #Should it be 0:(n*T - 1)%%(n) == i - 1, used to be 1:(n*T)%%n == i - 1
  A[i, which((0:(n * T - 1))%%n == i - 1)] <- 1
}



constr.st <- list(A = A, e = rep(0, dim(A)[1]))

#scaled_RW_prec <- inla.scale.model(struct_RW1,
#                                   list(A = matrix(1, 1, dim(struct_RW1)[1]),
#                                                    e = 0))

# Kronecker product between RW1 and IID space term
# order matters here! In our data set, the time ordering take precedent over the space ordering
# so we must have the RW on the left and space on the right
R <- struct_RW1 %x% diag(n)

interaction_param = c(1, 0.01)

typeII_formula <- update(basic_linear_formula,
                         ~. + f(space_time_unstructured, 
                                model = "generic0", 
                                Cmatrix = R, 
                                extraconstr = constr.st, 
                                rankdef = n, 
                                param = interaction_param))



ptm <- Sys.time()
modII <- inla(typeII_formula,
              data = ohio_df,
              family = "poisson",
              E = pop_at_risk,
              control.compute = list(config = TRUE,
                                     cpo = TRUE,
                                     waic = TRUE))
time = Sys.time() - ptm

#Check that constraints actually sum to zero over each RW for every county

print(c("Number of constraints (should be 2 + 88 = 90): ",toString(modII$misc$configs$constr$nc)))
print(c("mean CPO:", toString(round(mean(modII$cpo$cpo), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(modII$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(modII$cpo$failure[1])))
print(c("CPU: ", toString(summary(modII)$cpu.used)))
print(c("Time: ", time))
###

#Fit model w. type III interaction
###
#Type III interaction

#- Constraints: The ICAR at each time needs to sum to 0

### creating constraint needed for Type III Interaction
A <- matrix(0, nrow = T, ncol = n * T)
for (i in 1:T) {
  # The ICAR at each time point needs to sum to 0
  A[i, ((i - 1) * n + 1):(i * n)] <- 1
}

### defining Kronecker product for the Type III Interaction
# get scaled ICAR
#scaled_ICAR_prec <- INLA::inla.scale.model(ICAR_structure, 
#                                           constr = list(A = matrix(1,
#                                                                    1,
#                                                                    dim(ICAR_structure)[1]),
#                                                         e = 0))

# Kronecker product between IID x ICAR
R <- diag(T) %x% ICAR_structure 

# specify constraints in INLA-ready format
constr.st <- list(A = A, e = rep(0, dim(A)[1]))

interactions_param = c(1, 0.01)

#Model w. type III interaction
typeIII_formula <- update(basic_linear_formula, 
                          ~. + f(space_time_unstructured, 
                                 model = "generic0", 
                                 Cmatrix = R, 
                                 extraconstr = constr.st, 
                                 rankdef = T, 
                                 param = interaction_param))


ptm <- Sys.time()
modIII <- inla(typeIII_formula,
               data = ohio_df,
               family = "poisson",
               E = pop_at_risk,
               control.compute = list(config = TRUE, 
                                      cpo = TRUE,
                                      waic = TRUE))
time = Sys.time() - ptm


print(c("Number of constraints (should be 2 + 21 = 23): ",toString(modIII$misc$configs$constr$nc)))
print(c("mean CPO:", toString(round(mean(modIII$cpo$cpo), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(modIII$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(modIII$cpo$failure[1])))
print(c("CPU: ", toString(summary(modIII)$cpu.used)))
print(c("Time: ", time))
###

#Fit model w. type IV interaction
###
# Type IV Interaction ----------------------------------------------------------

# Specifying constraints needed for Type IV Interaction
time_constr <- matrix(0, n, n * T)
for (i in 1:n) {
  time_constr[i, which((0:(n * T - 1))%%n == i - 1)] <- 1
}
space_constr <- matrix(0, T-1, n * T)
# theoretically it's enough to only go through N-1 here
# note that we could have instead done 1:(T-1) in the first loop and 1:n in the second
for (i in 1:(T-1)) { 
  space_constr[i, ((i - 1) * n + 1):(i * n)] <- 1
}


tmp <- rbind(time_constr, space_constr) 
constr.st <- list(A = tmp, e = rep(0, dim(tmp)[1]))

# Kronecker product between ICAR x RW1
R <- struct_RW1 %x% ICAR_structure

formulaIV <- update(basic_linear_formula,
                    ~. + f(space_time_unstructured, 
                           model = "generic0",
                           Cmatrix = R,
                           extraconstr = constr.st,
                           rankdef = n + T - 1, 
                           param = c(1, 0.01)))


print("Done diddily dooo")

#
#ptm <- Sys.time()
modIV <- inla(formulaIV,
              data = ohio_df,
              family = "poisson",
              E = pop_at_risk,
              control.compute = list(config = TRUE,
                                     cpo = TRUE,
                                     waic = TRUE))

print("Gurkaaaaaaaaaaaa")

#time = Sys.time() - ptm
###