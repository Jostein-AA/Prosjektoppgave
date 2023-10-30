################################################################################
#Libraries and environment preparation
#Clear environment
rm(list = ls())

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

#Load other files
if(getwd() != "C:/Users/joste/Documents/H2023/Code/Prosjektoppgave"){
  setwd("H2023/Code/Prosjektoppgave/")
}

source("functions_plotting_etc.R")
################################################################################
#Load in necessities
#Load structure matrices
struct_RW1 = readMM(file='struct_RW1.txt')
struct_RW2 = readMM(file = 'struct_RW2.txt')
ICAR_structure = readMM(file = 'ICAR_struct.txt')

#Load data
ohio_df <- read.csv("ohio_df.csv")

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]

#Get number of counties (n) and number of years (T)
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points
################################################################################
#Merge data with map for plotting over the areas
ohio_df_w_geometries <- merge(ohio_map, ohio_df,
                              by.x = c("NAME"), by.y = c("county_name"),
                              all = T, suffixes = T)


################################################################################
#Explorative plots

####### 
#Median rate for each county

#Start by finding median rate for each county
ohio_map$median_rate <- rep(0, n)  #Initialize to zero
for(county_name in ohio_df$county_name){ #Iterate over all the counties
  #Extract values for county 'county_name'
  temp <- ohio_df[ohio_df$county_name == county_name, ]
  
  #Calculate median rate for that county
  ohio_map$median_rate[ohio_map$NAME == county_name] = median(temp$rate)
}

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
#scale = heat.colors(8, rev= TRUE)
median_rate_plot <- ggplot(data = ohio_map) +  
  geom_sf(aes(fill = median_rate), #Plots cases per thousand
          alpha = 1,
          color="black") + ggtitle("Median rate for each county") + 
  theme(plot.title = element_text(size = 12),
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
median_rate_plot
#######
## Heatmap of rate over all years

#Calcualte min and max rate over all to help create one heatmap scale common for all
min_rate <- min(ohio_df$rate); max_rate <- max(ohio_df$rate)
hardcoded_bins = seq(min_rate, max_rate, length.out = 8)
hardcoded_bins = round(hardcoded_bins, 4)

p_1 <- case_count_plot_1_year(ohio_df_w_geometries, 1, hardcoded_bins = hardcoded_bins)
p_2 <- case_count_plot_1_year(ohio_df_w_geometries, 2, hardcoded_bins = hardcoded_bins)
p_3 <- case_count_plot_1_year(ohio_df_w_geometries, 3, hardcoded_bins = hardcoded_bins)
p_4 <- case_count_plot_1_year(ohio_df_w_geometries, 4, hardcoded_bins = hardcoded_bins)
p_5 <- case_count_plot_1_year(ohio_df_w_geometries, 5, hardcoded_bins = hardcoded_bins)
p_6 <- case_count_plot_1_year(ohio_df_w_geometries, 6, hardcoded_bins = hardcoded_bins)
p_7 <- case_count_plot_1_year(ohio_df_w_geometries, 7, hardcoded_bins = hardcoded_bins)
p_8 <- case_count_plot_1_year(ohio_df_w_geometries, 8, hardcoded_bins = hardcoded_bins)
p_9 <- case_count_plot_1_year(ohio_df_w_geometries, 9, hardcoded_bins = hardcoded_bins)
p_10 <- case_count_plot_1_year(ohio_df_w_geometries, 10, hardcoded_bins = hardcoded_bins)
p_11 <- case_count_plot_1_year(ohio_df_w_geometries, 11, hardcoded_bins = hardcoded_bins)
p_12 <- case_count_plot_1_year(ohio_df_w_geometries, 12, hardcoded_bins = hardcoded_bins)
p_13 <- case_count_plot_1_year(ohio_df_w_geometries, 13, hardcoded_bins = hardcoded_bins)
p_14 <- case_count_plot_1_year(ohio_df_w_geometries, 14, hardcoded_bins = hardcoded_bins)
p_15 <- case_count_plot_1_year(ohio_df_w_geometries, 15, hardcoded_bins = hardcoded_bins)
p_16 <- case_count_plot_1_year(ohio_df_w_geometries, 16, hardcoded_bins = hardcoded_bins)
p_17 <- case_count_plot_1_year(ohio_df_w_geometries, 17, hardcoded_bins = hardcoded_bins)
p_18 <- case_count_plot_1_year(ohio_df_w_geometries, 18, hardcoded_bins = hardcoded_bins)
p_19 <- case_count_plot_1_year(ohio_df_w_geometries, 19, hardcoded_bins = hardcoded_bins)
p_20 <- case_count_plot_1_year(ohio_df_w_geometries, 20, hardcoded_bins = hardcoded_bins)
p_21 <- case_count_plot_1_year(ohio_df_w_geometries, 21, hardcoded_bins = hardcoded_bins)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")


#######
#Depth measure over time series
#Transform data to functional data format (only for calculating depth measures)
grid = unique(ohio_df$year)
temp_ = ohio_df[ohio_df$county_name == unique(ohio_df$county_name)[1], ]
temp_data = temp_$rate
for(county in unique(ohio_df$county_name)[2:length(unique(ohio_df$county_name))]){
  temp_ = ohio_df[ohio_df$county_name == county, ]
  temp_data = rbind(temp_data, temp_$rate)
}
fdata <- fData(grid, temp_data)

#Calculate the median curve by way of modified band depth
median_curve = data.frame(
                 median_curve = median_fData(fdata, type = "MBD")$values[1, ],
                 year = grid)


#Plot the functional data (each county has its own curve)
#along with the depth median to see potential temporal trend in data
ggplot(data = ohio_df, 
       aes(x = year, 
           y = rate, 
           group = county_name, 
           color = county_name)) + 
  theme(legend.position="none") +
  geom_line() + 
  geom_line(data = median_curve, 
            aes(x = year, y = median_curve), 
            col = 'black', lwd = 3)

################################################################################
#Specify base linear formula and hyperparameters and corresponding priors
#for the temporal and spatial effects
######################################
#Specify hyperparameters with corresponding priors
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


######################################
#Make the base formula (is the overall mean specified correctly here???)
basic_linear_formula <- deaths ~ 1 + f(year, 
                                       model = 'bym',
                                       scale.model = T, 
                                       constr = T, 
                                       rankdef = 1,
                                       graph = struct_RW1,
                                       hyper = temporal_hyper) + 
                                     f(county, 
                                       model = 'bym',
                                       scale.model = T,
                                       constr = T,
                                       rankdef = 1,
                                       graph = ICAR_structure,
                                       hyper = spatial_hyper
                                       )
################################################################################
#Fit the basic model
######################################
#Fit the basic model w.o. interactions
#Measure time to do inference
ptm <- Sys.time()
basic_model_fit <- inla(basic_linear_formula,
                        data = ohio_df,
                        family = "poisson",
                        E = pop_at_risk, #Is it supposed to be expected or population at risk???
                        control.compute = list(config = TRUE, # needed if you want to see constraints later
                                               cpo = T,
                                               waic = T), #Needed for model selection
                        verbose = F #Not T unless want to see a lot of stuff
) 
time = Sys.time()-ptm


print(c("Number of constraints (should be 2): ",toString(basic_model_fit$misc$configs$constr$nc)))
print(c("- sum(log(CPO)):", toString(round(-sum(log(basic_model_fit$cpo$cpo)), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(basic_model_fit$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(basic_model_fit$cpo$failure[1])))
print(c("CPU: ", toString(summary(basic_model_fit)$cpu.used)))
print(c("Time to compute", toString(time)))

#When constr = T for a model, then the effects ought to sum-to-zero right?
#check it???

#Check sum-to-zero for temporal structured effect
print("Temporal structured effects sum: ")
print(sum(basic_model_fit$summary.random$year$mode[1:T]))
print(sum(basic_model_fit$summary.random$year$mode[(T+1):(2*T)]))
print(sum(basic_model_fit$summary.random$year$mode))


#Check sum-to-zero for spatial structured effect
print("Spatial structured effects sum: ")
print(sum(basic_model_fit$summary.random$county$mode[1:n]))
print(sum(basic_model_fit$summary.random$county$mode[(n+1):(2*n)]))
print(sum(basic_model_fit$summary.random$county$mode))

#Save cpo summary in easy way
base_cpo_summary <-  round(-sum(log(basic_model_fit$cpo$cpo)), digits = 4)

################
#Inference on the basic_model_fit
plot(basic_model_fit)

##########
#Plot the temporal effects

#Format results 
temporal_effects_fitted <- data.frame(year = 1:T)
temporal_effects_fitted$q025_struct_effect <- basic_model_fit$summary.random$year$'0.025quant'[1:T]
temporal_effects_fitted$median_struct_effect <- basic_model_fit$summary.random$year$'0.5quant'[1:T]
temporal_effects_fitted$q975_struct_effect <- basic_model_fit$summary.random$year$'0.975quant'[1:T]



#
ggplot(data = temporal_effects_fitted) + ggtitle("Median structured Temporal effect") +
  ylab("Temporal effect") + 
  geom_line(aes(x = year, y = q025_struct_effect), col = "red") + 
  geom_line(aes(x = year, y = q975_struct_effect), col = "blue") +
  geom_ribbon(aes(x = year, ymin = q025_struct_effect, ymax = q975_struct_effect),
              fill = "grey70") + 
  geom_line(aes(x = year, y = median_struct_effect)) + 
  xlim(1, 21)



##########
#Plot the spatial effects
spatial_structured_effect_median <- basic_model_fit$summary.random$county$'0.5quant'[1:n]
spatial_structured_effect_q025 <- basic_model_fit$summary.random$county$'0.025quant'[1:n]
spatial_structured_effect_q975 <- basic_model_fit$summary.random$county$'0.975quant'[1:n]
spatial_structured_effect_sd   <- basic_model_fit$summary.random$county$sd[1:n]

temp_ohio_map <- ohio_map[ ,c("geometry", "NAME")]
temp_ohio_map$median <- spatial_structured_effect_median
temp_ohio_map$q025 <- spatial_structured_effect_q025
temp_ohio_map$q975 <- spatial_structured_effect_q975
temp_ohio_map$sd   <- spatial_structured_effect_sd

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale_1 = scale_col[c(3,7,10,14,18,21,25,30)] #Select color scale to be more red

p_1 <- ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = median), 
          alpha = 1,
          color="black") + ggtitle("Median Spatial Structured Effect each County") +
  theme(plot.title = element_text(size = 12),
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
    palette = function(x) c(scale_1),
    labels = function(x){x},
    guide = "colorscale")

#Plot spatial structured effect along with median rate side by side
ggarrange(p_1, median_rate_plot,
          ncol = 2, nrow = 1,
          common.legend = FALSE)
#Plots exhibit similar structure

ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = sd), 
          alpha = 1,
          color="black") + ggtitle("Standard Deviation of Spatial structured\n effect for each county") +
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

##########
#Plot the fitted values

#####
###Heatmap of fitted values

##Just the fitted values
#Extract fitted values
fitted_values <- data.frame(basic_model_fit$summary.fitted.values$mean)

#fitted values which is ???
colnames(fitted_values) <- "rate"
fitted_values$year <- ohio_df$year
fitted_values$county <- ohio_df$county
fitted_values$county_name <- ohio_df$county_name
fitted_values$to_sort_on <- ohio_df$space_time_unstructured

#Merge with geometry, this messes up order
fitted_values <- merge(ohio_map, fitted_values,
                       by.x = c("NAME"), by.y = c("county_name"),
                       all = T, suffixes = T)
#Sort on what is space_time_unstructured in ohio_df to get right order
fitted_values <- fitted_values[order(fitted_values$to_sort_on), ]

#Fix indices
rownames(fitted_values) <- 1:nrow(fitted_values)    # Assign sequence to row names

#Create heatmaps of the fitted values
p_1 <- case_count_plot_1_year(fitted_values, 1)
p_2 <- case_count_plot_1_year(fitted_values, 2)
p_3 <- case_count_plot_1_year(fitted_values, 3)
p_4 <- case_count_plot_1_year(fitted_values, 4)
p_5 <- case_count_plot_1_year(fitted_values, 5)
p_6 <- case_count_plot_1_year(fitted_values, 6)
p_7 <- case_count_plot_1_year(fitted_values, 7)
p_8 <- case_count_plot_1_year(fitted_values, 8)
p_9 <- case_count_plot_1_year(fitted_values, 9)
p_10 <- case_count_plot_1_year(fitted_values, 10)
p_11 <- case_count_plot_1_year(fitted_values, 11)
p_12 <- case_count_plot_1_year(fitted_values, 12)
p_13 <- case_count_plot_1_year(fitted_values, 13)
p_14 <- case_count_plot_1_year(fitted_values, 14)
p_15 <- case_count_plot_1_year(fitted_values, 15)
p_16 <- case_count_plot_1_year(fitted_values, 16)
p_17 <- case_count_plot_1_year(fitted_values, 17)
p_18 <- case_count_plot_1_year(fitted_values, 18)
p_19 <- case_count_plot_1_year(fitted_values, 19)
p_20 <- case_count_plot_1_year(fitted_values, 20)
p_21 <- case_count_plot_1_year(fitted_values, 21)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")

#The difference between the fitted values and actual values
#Make a copy to get same structure
fitted_values_copy = fitted_values
fitted_values_copy$rate = fitted_values_copy$rate - ohio_df$rate #abs()

min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")
#Output shows that for the most part, most fitted values for the rate is
#Quite close to the actual value


#Relative deviation
fitted_values_copy$rate = ifelse(ohio_df$rate == 0, 
                                 (fitted_values$rate)/0.0001,
                                 abs((fitted_values$rate - ohio_df$rate)/ohio_df$rate))
  


min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
#hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)
hardcoded_bins_diff = seq(0, 2, length.out = 8)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")


#####
##Timeseries of fitted values
#Just the fitted values

#The difference between the fitted values and actual values

################################################################################
#Fit the model w. type I interaction

####################################
# Type I Interaction hyperparameters and hyper prior
typeI_hyperparameters_priors = list(theta=list(prior="pc.prec",
                                               param=c(1,0.01)))

#Make the formula for the model w. type I interaction by updating base formula
typeI_formula <- update(basic_linear_formula, 
                        ~. + f(space_time_unstructured,
                               model="iid", #Has to be iid, whole point
                               hyper = typeI_hyperparameters_priors ))





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
print(c("-sum(log(CPO)): ", toString(round(-sum(log(typeI_fit$cpo$cpo)), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(typeI_fit$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(typeI_fit$cpo$failure[1])))
print(c("CPU: ", toString(summary(typeI_fit)$cpu.used)))
print(c("Time: ", time))

#Save type I cpo, notice that it is smaller than base (that is good)
typeI_cpo_summary <- round(-sum(log(typeI_fit$cpo$cpo)), digits = 4)

##################
#Inference on results

plot(typeI_fit)

######
#Plot the temporal effects

#Format results 
temporal_effects_fitted <- data.frame(year = 1:T)
temporal_effects_fitted$q025_struct_effect <- typeI_fit$summary.random$year$'0.025quant'[1:T]
temporal_effects_fitted$median_struct_effect <- typeI_fit$summary.random$year$'0.5quant'[1:T]
temporal_effects_fitted$q975_struct_effect <- typeI_fit$summary.random$year$'0.975quant'[1:T]



#
ggplot(data = temporal_effects_fitted) + ggtitle("Median structured Temporal effect") +
  ylab("Temporal effect") + 
  geom_line(aes(x = year, y = q025_struct_effect), col = "red") + 
  geom_line(aes(x = year, y = q975_struct_effect), col = "blue") +
  geom_ribbon(aes(x = year, ymin = q025_struct_effect, ymax = q975_struct_effect),
              fill = "grey70") + 
  geom_line(aes(x = year, y = median_struct_effect)) + 
  xlim(1, 21)


######
#Plot the spatial effects
spatial_structured_effect_median <- typeI_fit$summary.random$county$'0.5quant'[1:n]
spatial_structured_effect_q025 <- typeI_fit$summary.random$county$'0.025quant'[1:n]
spatial_structured_effect_q975 <- typeI_fit$summary.random$county$'0.975quant'[1:n]
spatial_structured_effect_sd   <- typeI_fit$summary.random$county$sd[1:n]

temp_ohio_map <- ohio_map[ ,c("geometry", "NAME")]
temp_ohio_map$median <- spatial_structured_effect_median
temp_ohio_map$q025 <- spatial_structured_effect_q025
temp_ohio_map$q975 <- spatial_structured_effect_q975
temp_ohio_map$sd   <- spatial_structured_effect_sd

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale_1 = scale_col[c(3,7,10,14,18,21,25,30)] #Select color scale to be more red

p_1 <- ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = median), 
          alpha = 1,
          color="black") + ggtitle("Median Spatial Structured Effect each County") +
  theme(plot.title = element_text(size = 12),
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
    palette = function(x) c(scale_1),
    labels = function(x){x},
    guide = "colorscale")

#Plot median structured spatial effect in each county
#Side by side with median rate in each county
ggarrange(p_1, median_rate_plot,
          ncol = 2, nrow = 1,
          common.legend = FALSE)
#The plots exhibit a fairly similar structure, something which makes a lot of sense

ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = sd), 
          alpha = 1,
          color="black") + ggtitle("Standard Deviation of Spatial structured\n effect for each county") +
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

#######
#Plot the fitted values

#####
###Heatmap of fitted values

##Just the fitted values
#Extract fitted values
fitted_values <- data.frame(typeI_fit$summary.fitted.values$mean)

#fitted values which is ???
colnames(fitted_values) <- "rate"
fitted_values$year <- ohio_df$year
fitted_values$county <- ohio_df$county
fitted_values$county_name <- ohio_df$county_name
fitted_values$to_sort_on <- ohio_df$space_time_unstructured

#Merge with geometry, this messes up order
fitted_values <- merge(ohio_map, fitted_values,
                       by.x = c("NAME"), by.y = c("county_name"),
                       all = T, suffixes = T)
#Sort on what is space_time_unstructured in ohio_df to get right order
fitted_values <- fitted_values[order(fitted_values$to_sort_on), ]

#Fix indices
rownames(fitted_values) <- 1:nrow(fitted_values)    # Assign sequence to row names

#Create heatmaps of the fitted values
p_1 <- case_count_plot_1_year(fitted_values, 1)
p_2 <- case_count_plot_1_year(fitted_values, 2)
p_3 <- case_count_plot_1_year(fitted_values, 3)
p_4 <- case_count_plot_1_year(fitted_values, 4)
p_5 <- case_count_plot_1_year(fitted_values, 5)
p_6 <- case_count_plot_1_year(fitted_values, 6)
p_7 <- case_count_plot_1_year(fitted_values, 7)
p_8 <- case_count_plot_1_year(fitted_values, 8)
p_9 <- case_count_plot_1_year(fitted_values, 9)
p_10 <- case_count_plot_1_year(fitted_values, 10)
p_11 <- case_count_plot_1_year(fitted_values, 11)
p_12 <- case_count_plot_1_year(fitted_values, 12)
p_13 <- case_count_plot_1_year(fitted_values, 13)
p_14 <- case_count_plot_1_year(fitted_values, 14)
p_15 <- case_count_plot_1_year(fitted_values, 15)
p_16 <- case_count_plot_1_year(fitted_values, 16)
p_17 <- case_count_plot_1_year(fitted_values, 17)
p_18 <- case_count_plot_1_year(fitted_values, 18)
p_19 <- case_count_plot_1_year(fitted_values, 19)
p_20 <- case_count_plot_1_year(fitted_values, 20)
p_21 <- case_count_plot_1_year(fitted_values, 21)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")

#From the plot it seems like maybe the fitted values go in a 'circular' fashion
#Rates at time 1 seem similar to that at time 21
#May be due to non hardcoded bins...


#The difference between the fitted values and actual values
#Make a copy to get same structure
fitted_values_copy = fitted_values
fitted_values_copy$rate = fitted_values_copy$rate - ohio_df$rate #abs()

min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")
#Output shows that for the most part, most fitted values for the rate is
#Quite close to the actual value


#Relative deviation
fitted_values_copy$rate = ifelse(ohio_df$rate == 0, 
                                 (fitted_values$rate)/0.0001,
                                 abs((fitted_values$rate - ohio_df$rate)/ohio_df$rate))



min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
#hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)
hardcoded_bins_diff = seq(0, 2, length.out = 8)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")
#Looks pretty good actually, some regions at early times that are quite far off, 
#But other than that we decently good

################################################################################
#Type II

##################################################
# creating constraint matrix needed for Type II Interaction

#sum-to-zero constraint be wilding (Is it the right way??? Does it actually do what it is supposed to???)
# - n rows as this interaction makes it so there is a RW for each county independent of each other
# - n * T columns as there are a total of n * T interaction terms (delta_(it))
# - sets A[i, which((1:(n * T))%%n == i)] = 1, as for county i, the interaction (delta_(it))
#     is only dependent on (delta_(i,t+-1)). Therefore the sum-to-zero is over these 88 RWs
#     additionaly, set A[n, which(1:(n * T)%%n == 0)] <- 1 
# - Constraints: The RW1 in each area needs to sum to 0. Hence in constr.st e is a zero vector
A <- matrix(0, nrow = n, ncol = n * T)
for (i in 1:(n - 1)) {
  #I disagree with this constraint
  #A[i, which((1:(n * T))%%n == i - 1)] <- 1
  
  #How do I actually check that this works?
  A[i, which((1:(n * T))%%n == i)] <- 1
}
A[n, which((1:(n * T))%%n == 0)] <- 1


constr.st <- list(A = A, e = rep(0, dim(A)[1]))

scaled_RW_prec <- inla.scale.model(struct_RW1,
                                   list(A = matrix(1, 1, dim(struct_RW1)[1]),
                                                    e = 0))

##################################################
#Fit the model w. type II interaction

# Kronecker product between RW1 and IID space term
# order matters here! In our data set, the time ordering take precedent over the space ordering
# so we must have the RW on the left and space on the right
R <- scaled_RW_prec %x% diag(n)

typeII_hyperparameters_priors = list(theta=list(prior="pc.prec",
                                                param=c(1,0.01)))

typeII_formula <- update(basic_linear_formula,
                         ~. + f(space_time_unstructured, 
                                model = "generic0", 
                                Cmatrix = R, 
                                extraconstr = constr.st, 
                                rankdef = n, 
                                hyper = typeII_hyperparameters_priors))



ptm <- Sys.time()
typeII_fit <- inla(typeII_formula,
                   data = ohio_df,
                   family = "poisson",
                   E = pop_at_risk,
                   control.compute = list(config = TRUE,
                                          cpo = TRUE,
                                          waic = TRUE),
                   control.predictor = list(compute = TRUE)
                   )
time = Sys.time() - ptm

print(c("Number of constraints (should be 2 + 88 = 90): ",toString(typeII_fit$misc$configs$constr$nc)))
print(c("-sum(log(CPO)): ", toString(round(-sum(log(typeII_fit$cpo$cpo)), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(typeII_fit$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(typeII_fit$cpo$failure[1])))
print(c("CPU: ", toString(summary(typeII_fit)$cpu.used)))
print(c("Time: ", time))


#Check sum-to-zero constraints on interactions
interaction_mean <- typeII_fit$summary.random$space_time_unstructured$mean

print(constr.st$A %*% interaction_mean)

check_sum_to_zero <- ifelse(constr.st$A %*% interaction_mean > 0 - 1E-10 & 
                            constr.st$A %*% interaction_mean < 0 + 1E-10, 
                            0,
                            1)

print(c("Number of counties for which sum-to-zero did not hold: ", toString(sum(check_sum_to_zero))))


#save type II cpo
typeII_cpo_summary <- round(-sum(log(typeII_fit$cpo$cpo)), digits = 4)
#Notice that it is smaller than both typeI and hence also base (could be good)
#From looking at the time series for each county this makes some sense,
#as there are points in time where the rate varies a lot for certain counties

##################
#Inference on results

plot(typeII_fit)
#Why is interaction weird af???

######
#Plot the temporal effects

#Format results 
temporal_effects_fitted <- data.frame(year = 1:T)
temporal_effects_fitted$q025_struct_effect <- typeII_fit$summary.random$year$'0.025quant'[1:T]
temporal_effects_fitted$median_struct_effect <- typeII_fit$summary.random$year$'0.5quant'[1:T]
temporal_effects_fitted$q975_struct_effect <- typeII_fit$summary.random$year$'0.975quant'[1:T]



#
ggplot(data = temporal_effects_fitted) + ggtitle("Median structured Temporal effect") +
  ylab("Temporal effect") + 
  geom_line(aes(x = year, y = q025_struct_effect), col = "red") + 
  geom_line(aes(x = year, y = q975_struct_effect), col = "blue") +
  geom_ribbon(aes(x = year, ymin = q025_struct_effect, ymax = q975_struct_effect),
              fill = "grey70") + 
  geom_line(aes(x = year, y = median_struct_effect)) + 
  xlim(1, 21)


######
#Plot the spatial effects
spatial_structured_effect_median <- typeII_fit$summary.random$county$'0.5quant'[1:n]
spatial_structured_effect_q025 <- typeII_fit$summary.random$county$'0.025quant'[1:n]
spatial_structured_effect_q975 <- typeII_fit$summary.random$county$'0.975quant'[1:n]
spatial_structured_effect_sd   <- typeII_fit$summary.random$county$sd[1:n]

temp_ohio_map <- ohio_map[ ,c("geometry", "NAME")]
temp_ohio_map$median <- spatial_structured_effect_median
temp_ohio_map$q025 <- spatial_structured_effect_q025
temp_ohio_map$q975 <- spatial_structured_effect_q975
temp_ohio_map$sd   <- spatial_structured_effect_sd

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale_1 = scale_col[c(3,7,10,14,18,21,25,30)] #Select color scale to be more red

p_1 <- ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = median), 
          alpha = 1,
          color="black") + ggtitle("Median Spatial Structured Effect each County") +
  theme(plot.title = element_text(size = 12),
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
    palette = function(x) c(scale_1),
    labels = function(x){x},
    guide = "colorscale")

#Plot median structured spatial effect in each county
#Side by side with median rate in each county
ggarrange(p_1, median_rate_plot,
          ncol = 2, nrow = 1,
          common.legend = FALSE)
#The plots exhibit a fairly similar structure, something which makes a lot of sense

ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = sd), 
          alpha = 1,
          color="black") + ggtitle("Standard Deviation of Spatial structured\n effect for each county") +
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

#######
#Plot the fitted values

#####
###Heatmap of fitted values

##Just the fitted values
#Extract fitted values
fitted_values <- data.frame(typeII_fit$summary.fitted.values$mean)

#fitted values which is ???
colnames(fitted_values) <- "rate"
fitted_values$year <- ohio_df$year
fitted_values$county <- ohio_df$county
fitted_values$county_name <- ohio_df$county_name
fitted_values$to_sort_on <- ohio_df$space_time_unstructured

#Merge with geometry, this messes up order
fitted_values <- merge(ohio_map, fitted_values,
                       by.x = c("NAME"), by.y = c("county_name"),
                       all = T, suffixes = T)
#Sort on what is space_time_unstructured in ohio_df to get right order
fitted_values <- fitted_values[order(fitted_values$to_sort_on), ]

#Fix indices
rownames(fitted_values) <- 1:nrow(fitted_values)    # Assign sequence to row names

#Create heatmaps of the fitted values
p_1 <- case_count_plot_1_year(fitted_values, 1)
p_2 <- case_count_plot_1_year(fitted_values, 2)
p_3 <- case_count_plot_1_year(fitted_values, 3)
p_4 <- case_count_plot_1_year(fitted_values, 4)
p_5 <- case_count_plot_1_year(fitted_values, 5)
p_6 <- case_count_plot_1_year(fitted_values, 6)
p_7 <- case_count_plot_1_year(fitted_values, 7)
p_8 <- case_count_plot_1_year(fitted_values, 8)
p_9 <- case_count_plot_1_year(fitted_values, 9)
p_10 <- case_count_plot_1_year(fitted_values, 10)
p_11 <- case_count_plot_1_year(fitted_values, 11)
p_12 <- case_count_plot_1_year(fitted_values, 12)
p_13 <- case_count_plot_1_year(fitted_values, 13)
p_14 <- case_count_plot_1_year(fitted_values, 14)
p_15 <- case_count_plot_1_year(fitted_values, 15)
p_16 <- case_count_plot_1_year(fitted_values, 16)
p_17 <- case_count_plot_1_year(fitted_values, 17)
p_18 <- case_count_plot_1_year(fitted_values, 18)
p_19 <- case_count_plot_1_year(fitted_values, 19)
p_20 <- case_count_plot_1_year(fitted_values, 20)
p_21 <- case_count_plot_1_year(fitted_values, 21)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")

#From the plot it seems like maybe the fitted values go in a 'circular' fashion
#Rates at time 1 seem similar to that at time 21
#May be due to non hardcoded bins...


#The difference between the fitted values and actual values
#Make a copy to get same structure
fitted_values_copy = fitted_values
fitted_values_copy$rate = fitted_values_copy$rate - ohio_df$rate #abs()

min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")
#Output shows that for the most part, most fitted values for the rate is
#Quite close to the actual value


#Relative deviation
fitted_values_copy$rate = ifelse(ohio_df$rate == 0, 
                                 (fitted_values$rate)/0.0001,
                                 abs((fitted_values$rate - ohio_df$rate)/ohio_df$rate))



min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
#hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)
hardcoded_bins_diff = seq(0, 1, length.out = 8)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")
#Looks pretty good actually, some regions at early times that are quite far off, 
#But other than that we decently good

################################################################################
#Type III

##################################################
### creating constraint needed for Type III Interaction
A <- matrix(0, nrow = T, ncol = n * T)
for (i in 1:T) {
  # The ICAR at each time point needs to sum to 0
  A[i, ((i - 1) * n + 1):(i * n)] <- 1
}

# specify constraints in INLA-ready format
constr.st <- list(A = A, e = rep(0, dim(A)[1]))

### defining Kronecker product for the Type III Interaction
# get scaled ICAR
scaled_ICAR_prec <- INLA::inla.scale.model(ICAR_structure, 
                                           constr = list(A = matrix(1,
                                                                    1,
                                                                    dim(ICAR_structure)[1]),
                                                         e = 0))

# Kronecker product between IID x ICAR
R <- diag(T) %x% scaled_ICAR_prec 

##################################################
#Fit the model w. type III interaction


typeIII_hyperparameters_priors = list(theta=list(prior="pc.prec",
                                                param=c(1,0.01)))

#Model w. type III interaction
typeIII_formula <- update(basic_linear_formula, 
                          ~. + f(space_time_unstructured, 
                                 model = "generic0", 
                                 Cmatrix = R, 
                                 extraconstr = constr.st, 
                                 rankdef = T, 
                                 hyper = typeIII_hyperparameters_priors))



ptm <- Sys.time()
typeIII_fit <- inla(typeIII_formula,
                    data = ohio_df,
                    family = "poisson",
                    E = pop_at_risk,
                    control.compute = list(config = TRUE, 
                                           cpo = TRUE,
                                           waic = TRUE))
time = Sys.time() - ptm


print(c("Number of constraints (should be 2 + 21 = 23): ",toString(typeIII_fit$misc$configs$constr$nc)))
print(c("-sum(log(CPO)): ", toString(round(-sum(log(typeIII_fit$cpo$cpo)), digits = 4))))
print(c("CPO failiure, should have 1 unique value",
        toString(length(unique(typeIII_fit$cpo$failure))))) #-> soinla.cpo(res)
print(c("And that value is (should be 0): ", toString(typeIII_fit$cpo$failure[1])))
print(c("CPU: ", toString(summary(typeIII_fit)$cpu.used)))
print(c("Time: ", time))


#Check sum-to-zero constraints on interactions
interaction_mean <- typeIII_fit$summary.random$space_time_unstructured$mean

print(constr.st$A %*% interaction_mean)

check_sum_to_zero <- ifelse(constr.st$A %*% interaction_mean > 0 - 1E-10 & 
                              constr.st$A %*% interaction_mean < 0 + 1E-10, 
                            0,
                            1)

print(c("Number of counties for which sum-to-zero did not hold: ", toString(sum(check_sum_to_zero))))


#save type III cpo
typeIII_cpo_summary <- round(-sum(log(typeIII_fit$cpo$cpo)), digits = 4)

##################
#Inference on results

plot(typeIII_fit)
#Why is interaction weird af???

######
#Plot the temporal effects

#Format results 
temporal_effects_fitted <- data.frame(year = 1:T)
temporal_effects_fitted$q025_struct_effect <- typeIII_fit$summary.random$year$'0.025quant'[1:T]
temporal_effects_fitted$median_struct_effect <- typeIII_fit$summary.random$year$'0.5quant'[1:T]
temporal_effects_fitted$q975_struct_effect <- typeIII_fit$summary.random$year$'0.975quant'[1:T]



#
ggplot(data = temporal_effects_fitted) + ggtitle("Median structured Temporal effect") +
  ylab("Temporal effect") + 
  geom_line(aes(x = year, y = q025_struct_effect), col = "red") + 
  geom_line(aes(x = year, y = q975_struct_effect), col = "blue") +
  geom_ribbon(aes(x = year, ymin = q025_struct_effect, ymax = q975_struct_effect),
              fill = "grey70") + 
  geom_line(aes(x = year, y = median_struct_effect)) + 
  xlim(1, 21)


######
#Plot the spatial effects
spatial_structured_effect_median <- typeIII_fit$summary.random$county$'0.5quant'[1:n]
spatial_structured_effect_q025 <- typeIII_fit$summary.random$county$'0.025quant'[1:n]
spatial_structured_effect_q975 <- typeIII_fit$summary.random$county$'0.975quant'[1:n]
spatial_structured_effect_sd   <- typeIII_fit$summary.random$county$sd[1:n]

temp_ohio_map <- ohio_map[ ,c("geometry", "NAME")]
temp_ohio_map$median <- spatial_structured_effect_median
temp_ohio_map$q025 <- spatial_structured_effect_q025
temp_ohio_map$q975 <- spatial_structured_effect_q975
temp_ohio_map$sd   <- spatial_structured_effect_sd

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale_1 = scale_col[c(3,7,10,14,18,21,25,30)] #Select color scale to be more red

p_1 <- ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = median), 
          alpha = 1,
          color="black") + ggtitle("Median Spatial Structured Effect each County") +
  theme(plot.title = element_text(size = 12),
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
    palette = function(x) c(scale_1),
    labels = function(x){x},
    guide = "colorscale")

#Plot median structured spatial effect in each county
#Side by side with median rate in each county
ggarrange(p_1, median_rate_plot,
          ncol = 2, nrow = 1,
          common.legend = FALSE)
#The plots exhibit a fairly similar structure, something which makes a lot of sense

ggplot(data = temp_ohio_map) + 
  geom_sf(aes(fill = sd), 
          alpha = 1,
          color="black") + ggtitle("Standard Deviation of Spatial structured\n effect for each county") +
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

#######
#Plot the fitted values

#####
###Heatmap of fitted values

##Just the fitted values
#Extract fitted values
fitted_values <- data.frame(typeIII_fit$summary.fitted.values$mean)

#fitted values which is ???
colnames(fitted_values) <- "rate"
fitted_values$year <- ohio_df$year
fitted_values$county <- ohio_df$county
fitted_values$county_name <- ohio_df$county_name
fitted_values$to_sort_on <- ohio_df$space_time_unstructured

#Merge with geometry, this messes up order
fitted_values <- merge(ohio_map, fitted_values,
                       by.x = c("NAME"), by.y = c("county_name"),
                       all = T, suffixes = T)
#Sort on what is space_time_unstructured in ohio_df to get right order
fitted_values <- fitted_values[order(fitted_values$to_sort_on), ]

#Fix indices
rownames(fitted_values) <- 1:nrow(fitted_values)    # Assign sequence to row names

#Create heatmaps of the fitted values
p_1 <- case_count_plot_1_year(fitted_values, 1)
p_2 <- case_count_plot_1_year(fitted_values, 2)
p_3 <- case_count_plot_1_year(fitted_values, 3)
p_4 <- case_count_plot_1_year(fitted_values, 4)
p_5 <- case_count_plot_1_year(fitted_values, 5)
p_6 <- case_count_plot_1_year(fitted_values, 6)
p_7 <- case_count_plot_1_year(fitted_values, 7)
p_8 <- case_count_plot_1_year(fitted_values, 8)
p_9 <- case_count_plot_1_year(fitted_values, 9)
p_10 <- case_count_plot_1_year(fitted_values, 10)
p_11 <- case_count_plot_1_year(fitted_values, 11)
p_12 <- case_count_plot_1_year(fitted_values, 12)
p_13 <- case_count_plot_1_year(fitted_values, 13)
p_14 <- case_count_plot_1_year(fitted_values, 14)
p_15 <- case_count_plot_1_year(fitted_values, 15)
p_16 <- case_count_plot_1_year(fitted_values, 16)
p_17 <- case_count_plot_1_year(fitted_values, 17)
p_18 <- case_count_plot_1_year(fitted_values, 18)
p_19 <- case_count_plot_1_year(fitted_values, 19)
p_20 <- case_count_plot_1_year(fitted_values, 20)
p_21 <- case_count_plot_1_year(fitted_values, 21)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")

#From the plot it seems like maybe the fitted values go in a 'circular' fashion
#Rates at time 1 seem similar to that at time 21
#May be due to non hardcoded bins...


#The difference between the fitted values and actual values
#Make a copy to get same structure
fitted_values_copy = fitted_values
fitted_values_copy$rate = fitted_values_copy$rate - ohio_df$rate #abs()

min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")
#Output shows that for the most part, most fitted values for the rate is
#Quite close to the actual value


#Relative deviation
fitted_values_copy$rate = ifelse(ohio_df$rate == 0, 
                                 (fitted_values$rate)/0.0001,
                                 abs((fitted_values$rate - ohio_df$rate)/ohio_df$rate))



min_diff_rate = min(fitted_values_copy$rate)
max_diff_rate = max(fitted_values_copy$rate)
#hardcoded_bins_diff = round(seq(min_diff_rate, max_diff_rate, length.out = 8), 4)
hardcoded_bins_diff = seq(0, 1, length.out = 8)

p_1 <- case_count_plot_1_year(fitted_values_copy, 1, hardcoded_bins = hardcoded_bins_diff)
p_2 <- case_count_plot_1_year(fitted_values_copy, 2, hardcoded_bins = hardcoded_bins_diff)
p_3 <- case_count_plot_1_year(fitted_values_copy, 3, hardcoded_bins = hardcoded_bins_diff)
p_4 <- case_count_plot_1_year(fitted_values_copy, 4, hardcoded_bins = hardcoded_bins_diff)
p_5 <- case_count_plot_1_year(fitted_values_copy, 5, hardcoded_bins = hardcoded_bins_diff)
p_6 <- case_count_plot_1_year(fitted_values_copy, 6, hardcoded_bins = hardcoded_bins_diff)
p_7 <- case_count_plot_1_year(fitted_values_copy, 7, hardcoded_bins = hardcoded_bins_diff)
p_8 <- case_count_plot_1_year(fitted_values_copy, 8, hardcoded_bins = hardcoded_bins_diff)
p_9 <- case_count_plot_1_year(fitted_values_copy, 9, hardcoded_bins = hardcoded_bins_diff)
p_10 <- case_count_plot_1_year(fitted_values_copy, 10, hardcoded_bins = hardcoded_bins_diff)
p_11 <- case_count_plot_1_year(fitted_values_copy, 11, hardcoded_bins = hardcoded_bins_diff)
p_12 <- case_count_plot_1_year(fitted_values_copy, 12, hardcoded_bins = hardcoded_bins_diff)
p_13 <- case_count_plot_1_year(fitted_values_copy, 13, hardcoded_bins = hardcoded_bins_diff)
p_14 <- case_count_plot_1_year(fitted_values_copy, 14, hardcoded_bins = hardcoded_bins_diff)
p_15 <- case_count_plot_1_year(fitted_values_copy, 15, hardcoded_bins = hardcoded_bins_diff)
p_16 <- case_count_plot_1_year(fitted_values_copy, 16, hardcoded_bins = hardcoded_bins_diff)
p_17 <- case_count_plot_1_year(fitted_values_copy, 17, hardcoded_bins = hardcoded_bins_diff)
p_18 <- case_count_plot_1_year(fitted_values_copy, 18, hardcoded_bins = hardcoded_bins_diff)
p_19 <- case_count_plot_1_year(fitted_values_copy, 19, hardcoded_bins = hardcoded_bins_diff)
p_20 <- case_count_plot_1_year(fitted_values_copy, 20, hardcoded_bins = hardcoded_bins_diff)
p_21 <- case_count_plot_1_year(fitted_values_copy, 21, hardcoded_bins = hardcoded_bins_diff)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "right")


################################################################################
#Type IV

###################################################
#Make constraints









