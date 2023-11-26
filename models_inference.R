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

#Load in INLA objects
load("improper_RW1_ICAR_fitted.RData")
load("improper_RW2_ICAR_fitted.RData")




################################################################################
#Create latex table with results
#Format results data.frame for latex table
results.df <- data.frame(Model = rep(c("1", "2", "3", "4", "5"), 3),
                         RW = c(rep("1", 15), rep("2", 15)),
                         model_choice = rep(c("1", "2", "3"), 10),
                         value = 1:30)


#model = 1 -> base, model = 2 -> type I,..., model = 5 -> type IV
#RW = 1 -> RW1, RW = 2 -> RW2
#model_choice = 1 -> CPO, model_choice = 2 -> WAIC, model_choice = 3 -> Computational time (s)

#Insert values into table
#RW1 base model
results.df$value[results.df$RW == "1" & results.df$Model == "1" & results.df$model_choice == "1"] = mean(-log(RW1_ICAR_fit$cpo$cpo))
results.df$value[results.df$RW == "1" & results.df$Model == "1" & results.df$model_choice == "2"] = RW1_ICAR_fit$waic$waic
results.df$value[results.df$RW == "1" & results.df$Model == "1" & results.df$model_choice == "3"] = time_RW1_ICAR

#RW1 type I
results.df$value[results.df$RW == "1" & results.df$Model == "2" & results.df$model_choice == "1"] = mean(-log(RW1_ICAR_I_fit$cpo$cpo))
results.df$value[results.df$RW == "1" & results.df$Model == "2" & results.df$model_choice == "2"] = RW1_ICAR_I_fit$waic$waic
results.df$value[results.df$RW == "1" & results.df$Model == "2" & results.df$model_choice == "3"] = time_RW1_ICAR_I

#RW1 type II
results.df$value[results.df$RW == "1" & results.df$Model == "3" & results.df$model_choice == "1"] = mean(-log(RW1_ICAR_II_fit$cpo$cpo))
results.df$value[results.df$RW == "1" & results.df$Model == "3" & results.df$model_choice == "2"] = RW1_ICAR_II_fit$waic$waic
results.df$value[results.df$RW == "1" & results.df$Model == "3" & results.df$model_choice == "3"] = time_RW1_ICAR_II

#RW1 type III
results.df$value[results.df$RW == "1" & results.df$Model == "4" & results.df$model_choice == "1"] = mean(-log(RW1_ICAR_III_fit$cpo$cpo))
results.df$value[results.df$RW == "1" & results.df$Model == "4" & results.df$model_choice == "2"] = RW1_ICAR_III_fit$waic$waic
results.df$value[results.df$RW == "1" & results.df$Model == "4" & results.df$model_choice == "3"] = time_RW1_ICAR_III

#RW1 Type IV
results.df$value[results.df$RW == "1" & results.df$Model == "5" & results.df$model_choice == "1"] = mean(-log(RW1_ICAR_IV_fit$cpo$cpo))
results.df$value[results.df$RW == "1" & results.df$Model == "5" & results.df$model_choice == "2"] = RW1_ICAR_IV_fit$waic$waic
results.df$value[results.df$RW == "1" & results.df$Model == "5" & results.df$model_choice == "3"] = time_RW1_ICAR_IV

#RW2 base model
results.df$value[results.df$RW == "2" & results.df$Model == "1" & results.df$model_choice == "1"] = mean(-log(RW2_ICAR_fit$cpo$cpo))
results.df$value[results.df$RW == "2" & results.df$Model == "1" & results.df$model_choice == "2"] = RW2_ICAR_fit$waic$waic
results.df$value[results.df$RW == "2" & results.df$Model == "1" & results.df$model_choice == "3"] = time_RW2_ICAR

#RW2 type I
results.df$value[results.df$RW == "2" & results.df$Model == "2" & results.df$model_choice == "1"] = mean(-log(RW2_ICAR_I_fit$cpo$cpo))
results.df$value[results.df$RW == "2" & results.df$Model == "2" & results.df$model_choice == "2"] = RW2_ICAR_I_fit$waic$waic
results.df$value[results.df$RW == "2" & results.df$Model == "2" & results.df$model_choice == "3"] = time_RW2_ICAR_I

#RW2 type II
results.df$value[results.df$RW == "2" & results.df$Model == "3" & results.df$model_choice == "1"] = mean(-log(RW2_ICAR_II_fit$cpo$cpo))
results.df$value[results.df$RW == "2" & results.df$Model == "3" & results.df$model_choice == "2"] = RW2_ICAR_II_fit$waic$waic
results.df$value[results.df$RW == "2" & results.df$Model == "3" & results.df$model_choice == "3"] = time_RW2_ICAR_II

#RW2 type III
results.df$value[results.df$RW == "2" & results.df$Model == "4" & results.df$model_choice == "1"] = mean(-log(RW2_ICAR_III_fit$cpo$cpo))
results.df$value[results.df$RW == "2" & results.df$Model == "4" & results.df$model_choice == "2"] = RW2_ICAR_III_fit$waic$waic
results.df$value[results.df$RW == "2" & results.df$Model == "4" & results.df$model_choice == "3"] = time_RW2_ICAR_III

#RW2 Type IV
results.df$value[results.df$RW == "2" & results.df$Model == "5" & results.df$model_choice == "1"] = mean(-log(RW2_ICAR_IV_fit$cpo$cpo))
results.df$value[results.df$RW == "2" & results.df$Model == "5" & results.df$model_choice == "2"] = RW2_ICAR_IV_fit$waic$waic
results.df$value[results.df$RW == "2" & results.df$Model == "5" & results.df$model_choice == "3"] = time_RW2_ICAR_IV

#Make caption and label for latex table
caption = "Results"
label = "tab:all-results"

#Make latex table
latex_tabular <- latexTable(tabular(
  Heading()*RowFactor(RW, levelnames = c("RW1", "RW2"))*
    Heading("Model")*RowFactor(Model, levelnames = c("base model",
                                              "interaction type I",
                                              "interaction type II",
                                              "interaction type III",
                                              "interaction type IV"),
                        nopagebreak = "\\hline",
                        spacing = 0)~
    Heading()*Factor(model_choice, 
                     levelnames = c("CPO", "WAIC", "Computational time (s)"))*
    Heading()*value*Heading()*identity,
  data = results.df),
  caption = caption,
  label = label
)
latex_tabular

#Save latex table
cat(latex_tabular, file = "table.tex")

################################################################################


#Get the intercept and temporal structured effect plotted
plot_random_effects(RW1_ICAR_fit, ohio_map, n, T)

#Plot the posterior hyperparameters
plot_temporal_spatial_hyperparameters(RW1_ICAR_fit)











#####
#Code prison

#Use mean of cpo instead...
print_cpo_etc(RW1_ICAR_fit, time_RW1_ICAR)

#See a lot of things (intercept, precisions, effects, etc)
plot(RW1_ICAR_fit)

#Plot the intercept
plot_intercept(basic_model_fit)

#Plot posterior distributions of precision of random effects
plot_precisions_random_effects(basic_model_fit)

#See Structured temporal effect
plot_temporal_effect(basic_model_fit, T)

#See structured spatial effect plotted as heatmap
plot_spatial_effect(ohio_map, basic_model_fit, n)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(basic_model_fit, n, T)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, basic_model_fit)

###


#Inference on model w. type I interaction
print_cpo_etc(typeI_fit, time_typeI)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeI_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeI_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeI_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeI_fit, T)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeI_fit, n)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeI_fit, n, T)

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


#Use posterior samples to plot different ratios...



#Inference on type II model
print_cpo_etc(typeII_fit, time_typeII)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeII_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeII_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeII_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeII_fit, T)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeII_fit, n)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeII_fit, n, T)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeII_fit)


#Precision of interactions
plot_prec_interactions(typeII_fit, 
                       "Posterior density of precision of type II interaction")

#Interactions themselves
plot_interaction(typeII_fit)


###

#Type III

#Inference on type III model
print_cpo_etc(typeIII_fit, time_typeIII)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeIII_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeIII_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeIII_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeIII_fit, T)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeIII_fit, n)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeIII_fit, n, T)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeIII_fit)


#Precision of interactions
plot_prec_interactions(typeIII_fit, 
                       "Posterior density of precision of type III interaction")

#Interactions themselves
plot_interaction(typeIII_fit)

###

#Type IV

#Inference on type IV model
print_cpo_etc(typeIV_fit, time_typeIV)

#See a lot of things (intercept, precisions, effects, etc)
plot(typeIV_fit)

#Plot the intercept: Looks exactly like the one plotted for basic_model_fit (good)
plot_intercept(typeIV_fit)

#Plot posterior distributions of precision of random effects: Similar to basic_model_fit
plot_precisions_random_effects(typeIV_fit)

#See Structured temporal effect: Similar to basic_model_fit
plot_temporal_effect(typeIV_fit, T)

#See structured spatial effect plotted as heatmap: Similar to basic_model_fit
plot_spatial_effect(ohio_map, typeIV_fit, n)

#See fitted values for each county w. median, and 0.025%- and 0.975% quantiles
#along with actual values of rate as points
every_county_time_series(typeIV_fit, n, T)

#Plot the fitted values against the actual observed values
plot_fitted_vs_actual_together(ohio_df, typeIV_fit)


#Precision of interactions
plot_prec_interactions(typeIV_fit, 
                       "Posterior density of precision of type IV interaction")

#Interactions themselves
plot_interaction(typeIV_fit)

