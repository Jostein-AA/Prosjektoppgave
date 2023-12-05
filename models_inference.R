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
library(latex2exp)
library(tables)

#Set working directory
if(getwd() != "C:/Users/joste/Documents/H2023/Code/Prosjektoppgave"){
  setwd("H2023/Code/Prosjektoppgave/")
}
source("utilities.R")

#Load in INLA objects
load("improper_RW1_ICAR_fitted.RData"); load("improper_RW2_ICAR_fitted.RData")
load("proper_fitted.RData"); load("one_step_predictions.RData")



################################################################################
#Create latex table with results
#Format results data.frame for latex table
improper_results.df <- data.frame(Interaction = rep(c("1", "2", "3", "4", "5"), 3),
                                  RW = c(rep("1", 15), rep("2", 15)),
                                  model_choice = rep(c("1", "2", "3"), 10),
                                  value = 1:30)


#model = 1 -> base, model = 2 -> type I,..., model = 5 -> type IV
#RW = 1 -> RW1, RW = 2 -> RW2
#model_choice = 1 -> CPO, model_choice = 2 -> WAIC, model_choice = 3 -> Computational time (s)

#Insert values into table
#RW1 base model
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "1" & improper_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "1" & improper_results.df$model_choice == "2"] = RW1_ICAR_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "1" & improper_results.df$model_choice == "3"] = time_RW1_ICAR

#RW1 type I
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "2" & improper_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_I_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "2" & improper_results.df$model_choice == "2"] = RW1_ICAR_I_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "2" & improper_results.df$model_choice == "3"] = time_RW1_ICAR_I

#RW1 type II
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "3" & improper_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_II_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "3" & improper_results.df$model_choice == "2"] = RW1_ICAR_II_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "3" & improper_results.df$model_choice == "3"] = time_RW1_ICAR_II

#RW1 type III
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "4" & improper_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_III_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "4" & improper_results.df$model_choice == "2"] = RW1_ICAR_III_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "4" & improper_results.df$model_choice == "3"] = time_RW1_ICAR_III

#RW1 Type IV
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "5" & improper_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_IV_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "5" & improper_results.df$model_choice == "2"] = RW1_ICAR_IV_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "1" & improper_results.df$Interaction == "5" & improper_results.df$model_choice == "3"] = time_RW1_ICAR_IV

#RW2 base model
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "1" & improper_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "1" & improper_results.df$model_choice == "2"] = RW2_ICAR_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "1" & improper_results.df$model_choice == "3"] = time_RW2_ICAR

#RW2 type I
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "2" & improper_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_I_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "2" & improper_results.df$model_choice == "2"] = RW2_ICAR_I_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "2" & improper_results.df$model_choice == "3"] = time_RW2_ICAR_I

#RW2 type II
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "3" & improper_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_II_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "3" & improper_results.df$model_choice == "2"] = RW2_ICAR_II_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "3" & improper_results.df$model_choice == "3"] = time_RW2_ICAR_II

#RW2 type III
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "4" & improper_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_III_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "4" & improper_results.df$model_choice == "2"] = RW2_ICAR_III_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "4" & improper_results.df$model_choice == "3"] = time_RW2_ICAR_III

#RW2 Type IV
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "5" & improper_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_IV_fit$cpo$cpo))
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "5" & improper_results.df$model_choice == "2"] = RW2_ICAR_IV_fit$waic$waic
improper_results.df$value[improper_results.df$RW == "2" & improper_results.df$Interaction == "5" & improper_results.df$model_choice == "3"] = time_RW2_ICAR_IV

#Make caption and label for latex table
caption = "Reported values of CPO, WAIC, and computational time,
           for improper model specifications.
          First column indicates whether the structured temporal effect follows a RW1 or RW2.
          Second column states what type of interaction the model uses.
          Third column reports the calculated summary statistic of CPO.
          Fourth column gives the calculated WAIC for each model.
          Lastly the computational time is reported in the last column."
label = "improper-results"

#Make latex table
latex_tabular <- latexTable(tabular(
  Heading()*RowFactor(RW, levelnames = c("RW1", "RW2"))*
    Heading("Interaction")*RowFactor(Interaction, levelnames = c("No interaction",
                                              "Type I",
                                              "Type II",
                                              "Type III",
                                              "Type IV"),
                        nopagebreak = "\\hline",
                        spacing = 0)~
    Heading()*Factor(model_choice, 
                     levelnames = c("CPO", "WAIC", "Computational time (s)"))*
    Heading()*value*Heading()*identity,
  data = improper_results.df),
  caption = caption,
  label = label
)
latex_tabular

#Save latex table
cat(latex_tabular, file = "improper_table.tex")

#Make proper table
proper_results.df <- data.frame(Model = c(rep("1", 3), rep("2", 3), rep("3", 3)),
                                model_choice = rep(c("1", "2", "3"), 3),
                                value = 1:9)

# Model = 1 -> proper base model, Model = 2 -> proper interaction only, Model = 3 -> full proper model
#Model_choice = 1 -> CPO, 2 -> WAIC, 3 -> Computational time

#Proper base results
proper_results.df$value[proper_results.df$Model == "1" & proper_results.df$model_choice == "1"] = mean(-log(proper_base_fit$cpo$cpo))
proper_results.df$value[proper_results.df$Model == "1" & proper_results.df$model_choice == "2"] = proper_base_fit$waic$waic
proper_results.df$value[proper_results.df$Model == "1" & proper_results.df$model_choice == "3"] = time_proper_base

#Proper interaction results
proper_results.df$value[proper_results.df$Model == "2" & proper_results.df$model_choice == "1"] = mean(-log(proper_interaction_fit$cpo$cpo))
proper_results.df$value[proper_results.df$Model == "2" & proper_results.df$model_choice == "2"] = proper_interaction_fit$waic$waic
proper_results.df$value[proper_results.df$Model == "2" & proper_results.df$model_choice == "3"] = time_proper_interaction

#Full proper results
proper_results.df$value[proper_results.df$Model == "3" & proper_results.df$model_choice == "1"] = mean(-log(proper_full_fit$cpo$cpo))
proper_results.df$value[proper_results.df$Model == "3" & proper_results.df$model_choice == "2"] = proper_full_fit$waic$waic
proper_results.df$value[proper_results.df$Model == "3" & proper_results.df$model_choice == "3"] = time_proper_full * 60

#Make caption and label for latex table
caption_proper = "Reported values of CPO, WAIC, and computational time,
           for proper model specifications.
          First column states the model specification.
          Second column reports the calculated summary statistic of CPO.
          Third column gives the calculated WAIC for each model.
          Lastly the computational time is reported in the last column."
label_proper = "proper-results"


#Make table for proper results

#Make latex table
latex_tabular_proper <- latexTable(tabular(
    Heading("Model")*RowFactor(Model, levelnames = c("No interaction",
                                                     "Only interaction",
                                                     "Full model"),
                                     nopagebreak = "\\hline",
                                     spacing = 0)~
    Heading()*Factor(model_choice, 
                     levelnames = c("CPO", "WAIC", "Computational time (s)"))*
    Heading()*value*Heading()*identity,
  data = proper_results.df),
  caption = caption_proper,
  label = label_proper
)
latex_tabular_proper

#Save latex table
cat(latex_tabular_proper, file = "proper_table.tex")

################################################################################

############
##Plots of posteriors: save to pdf 7 by 3.2
plot_intercept(RW1_ICAR_fit, proper_base_fit)


#Plot the temporal random effect RW1 vs RW2: save to pdf 7.5 by 3.5
plot_temporal_effects_RW1_RW2(RW1_ICAR_fit, RW2_ICAR_fit, T)

#Plot temporal random effect ar1 + fixed (density fixed, ar1, ar1 + fixed)
#Save to pdf 9.5 by 3.5
plot_temporal_ar1(proper_base_fit)
#plot_temporal_ar1(proper_full_fit)


#Plot spatial effects Proper vs improper
plot_spatial_effects(RW1_ICAR_fit,
                     proper_base_fit, 
                     ohio_map,
                     n)


plot_spatial_std(RW1_ICAR_fit,
                 proper_base_fit, 
                 ohio_map,
                 n)




#Plot the posterior hyperparameters (should it include RW2?)
#Plot improper temporal posterior hyperparameters: save to pdf 7.5 by 3.5 
plot_improper_temporal_hyperparameters(RW1_ICAR_fit, RW2_ICAR_fit)


#Plot improper spatial posterior hyperparameters: 7 by 3.5
plot_improper_spatial_hyperparameters(RW1_ICAR_fit)


#Plot proper temporal hyperparameters: save to pdf 7.5 by 3.5
plot_proper_temporal_hyperparameter(proper_base_fit, proper_full_fit)


#Plot proper spatial hyperparameters: 7.5 by 3.5
plot_proper_spatial_hyperparameters(proper_base_fit, proper_full_fit)




#Plot the interactions
#improper: save to pdf 7.5 by 12.5
plot_improper_interaction(RW1_ICAR_I_fit,
                          RW1_ICAR_II_fit,
                          RW1_ICAR_III_fit,
                          RW1_ICAR_IV_fit,
                          RW2_ICAR_I_fit,
                          RW2_ICAR_II_fit,
                          RW2_ICAR_III_fit,
                          RW2_ICAR_IV_fit)


#Proper interaction
plot_proper_interaction(proper_interaction_fit, proper_full_fit)


#Plot hyperparameter for interactions
#RW1
plot_std_interactions_RW1(RW1_ICAR_I_fit,
                          RW1_ICAR_II_fit,
                          RW1_ICAR_III_fit,
                          RW1_ICAR_IV_fit)
#RW2
plot_std_interactions_RW2(RW2_ICAR_I_fit,
                          RW2_ICAR_II_fit,
                          RW2_ICAR_III_fit,
                          RW2_ICAR_IV_fit)


#Proper: save 10 by 3
plot_proper_hyperparameters(proper_interaction_fit,
                            proper_full_fit)

  
#####
#Plot fitted-values for certain regions against true values (as time series)

#Want to extract county with highest average rate, lowest average,
#the one with largest range, and one random one
max_avg_rate <- mean(ohio_df[ohio_df$county == 1, ]$rate)
id_max_avg_rate <- 1

min_avg_rate <- mean(ohio_df[ohio_df$county == 1, ]$rate)
id_min_avg_rate <- 1

max_range_rate <- max(ohio_df[ohio_df$county == 1, ]$rate) - min(ohio_df[ohio_df$county == 1, ]$rate)
id_max_range_rate <- 1
for(i in 2:n){
  county = ohio_df[ohio_df$county == i, ]
  curr_avg_rate <- mean(county$rate)
  curr_range <- max(county$rate) - min(county$rate)
  if(curr_avg_rate > max_avg_rate){
    max_avg_rate = curr_avg_rate
    id_max_avg_rate = i
  }
  if(curr_avg_rate < min_avg_rate){
    min_avg_rate = curr_avg_rate
    id_min_avg_rate = i
  }
  if(curr_range > max_range_rate){
    max_range_rate = curr_range
    id_max_range_rate = i
  }
}
print(id_max_avg_rate)
print(id_min_avg_rate)
print(id_max_range_rate)



counties = c(id_max_avg_rate, id_min_avg_rate, id_max_range_rate, 4)

#Produces plot of fitted vs true for
# first frame: county w. largest average rate
# second frame: county w. smallest average rate
# third frame: county w. largest range in rate
# fourth frame: a random county
select_county_timeseries(ohio_df,
                         RW1_ICAR_fit,
                         RW1_ICAR_II_fit,
                         proper_interaction_fit,
                         counties,
                         n,
                         T)


#somehow need to sort the fitted values for the proper models
proper_interaction_fit$summary.fitted.values

proper_interaction_fit$summary.random$county

#Just make sure that sorter of proper models work
#test <- sort_proper_fitted(ohio_df_changed, n, T)














#####
#Fitted values of best model against true values 
plot_fitted_vs_actual_together(ohio_df, RW1_ICAR_fit,
                               RW1_ICAR_II_fit, n, T)

#####
#Plot the heatmaps of some years (1968, 1973, 1978, 1983, 1988 maybe?) for some models
#Base model, RW1 type II, RW1 type IV, maybe proper ones
years_to_plot = c(1968, 1975, 1980, 1988)


#Extract true values for years_to_plot
if(ohio_df$year[1]==1){
  ohio_df$year = ohio_df$year + (1968 - 1)
}

actual <- ohio_df[ohio_df$year %in% years_to_plot, ]

#Merge true values with map
actual_n_map <- merge(ohio_map, actual,
                      by.x = c("NAME"), by.y = c("name"),
                      all = T, suffixes = T)

#Get death rate per 100 000
actual_n_map$rate <- actual_n_map$rate * 1E5 

#Hardcoded bins
hardcoded_bins = seq(min(actual_n_map$rate) - 0.5,
                     max(actual_n_map$rate) + 0.5,
                     length.out = 8)
hardcoded_bins = round(hardcoded_bins, 0)

p1 <- case_count_plot_1_year(actual_n_map, years_to_plot[1], hardcoded_bins)
p2 <- case_count_plot_1_year(actual_n_map, years_to_plot[2], hardcoded_bins)
p3 <- case_count_plot_1_year(actual_n_map, years_to_plot[3], hardcoded_bins)
p4 <- case_count_plot_1_year(actual_n_map, years_to_plot[4], hardcoded_bins)

ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1,
          common.legend = TRUE, legend = "right")


#Extract relative risk for base
#years_1_21 <- c(1, 8, 13, 21)
fitted_values_base_y1 <- RW1_ICAR_fit$summary.fitted.values$mean[1:(88)] * 1E5
fitted_values_base_y2 <- RW1_ICAR_fit$summary.fitted.values$mean[(7 * 88 + 1):(8 * 88)] * 1E5
fitted_values_base_y3 <- RW1_ICAR_fit$summary.fitted.values$mean[(12 * 88 + 1):(13 * 88)] * 1E5
fitted_values_base_y4 <- RW1_ICAR_fit$summary.fitted.values$mean[(20 * 88 + 1):(21 * 88)] * 1E5
min1 <- min(fitted_values_base_y1); min2 <- min(fitted_values_base_y2); min3<- min(fitted_values_base_y3); min4<- min(fitted_values_base_y4)
max1 <- max(fitted_values_base_y1); max2 <- max(fitted_values_base_y2); max3<- max(fitted_values_base_y3); max4<- max(fitted_values_base_y4)
#hardcoded_bins <- seq(min(min1, min2, min3, min4) - 1,
#                      max(max1, max2, max3, max4) + 1,
#                      length.out = 8)
#hardcoded_bins = round(hardcoded_bins, 0)

base_map <- ohio_map
base_map$rate <- fitted_values_base_y1; base_map$year = rep(1968, 88)
p1 <- case_count_plot_1_year(base_map, years_to_plot[1], hardcoded_bins)

base_map$rate <- fitted_values_base_y2; base_map$year = rep(1975, 88)
p2 <- case_count_plot_1_year(base_map, years_to_plot[2], hardcoded_bins)

base_map$rate <- fitted_values_base_y3; base_map$year = rep(1980, 88)
p3 <- case_count_plot_1_year(base_map, years_to_plot[3], hardcoded_bins)

base_map$rate <- fitted_values_base_y4; base_map$year = rep(1988, 88)
p4 <- case_count_plot_1_year(base_map, years_to_plot[4], hardcoded_bins)

ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1,
          common.legend = TRUE, legend = "right")

#####
#Violin plots of the relative risk
rates.df <- data.frame(fitted_rates = c(RW1_ICAR_fit$summary.fitted.values$mean,
                                        RW1_ICAR_I_fit$summary.fitted.values$mean,
                                        RW1_ICAR_II_fit$summary.fitted.values$mean,
                                        RW1_ICAR_III_fit$summary.fitted.values$mean,
                                        RW1_ICAR_IV_fit$summary.fitted.values$mean),
                       type = c(rep("base", length(RW1_ICAR_fit$summary.fitted.values$mean)),
                                rep("Type I", length(RW1_ICAR_I_fit$summary.fitted.values$mean)),
                                rep("Type II", length(RW1_ICAR_II_fit$summary.fitted.values$mean)),
                                rep("Type III", length(RW1_ICAR_III_fit$summary.fitted.values$mean)),
                                rep("Type IV", length(RW1_ICAR_IV_fit$summary.fitted.values$mean))))



violin_plot_rate(rates.df)


#Observed vs predicted plots as lines???



############
##One-step predictor ...

#Function to use inla.tmarginal in lapply

#my_inla_t_marginal <- function(prediction_marginal){
  #a function to use inla.tmarginal on several values at once
#  return(inla.tmarginal(function(x){exp(x)}, prediction_marginal))
#}

#base_predicted = lapply(base_prediction_marginals, FUN = my_inla_t_marginal)
#I_predicted = lapply(I_prediction_marginals, FUN = my_inla_t_marginal)
#II_predicted = lapply(II_prediction_marginals, FUN = my_inla_t_marginal)
#III_predicted = lapply(III_prediction_marginals, FUN = my_inla_t_marginal)
#IV_predicted = lapply(IV_prediction_marginals, FUN = my_inla_t_marginal)

crpsNormal <- function(x, mu = 0, sig = 1){
  ## Function to compute the CRPS under normality assumption
  ## Here: x denotes the actual observation and mu and sigma
  ## mean and sd of the predictive distribution.
  ## (see Held et al. (2010), page 1296
  
  x0 <- (x - mu) / sig
  res <- sig * (1 / sqrt(pi) -  2 * dnorm(x0) - x0 * (2 * pnorm(x0) - 1))
  
  ## sign as in Held (2008)
  res <- -res
  return(res)
}



#Do NOT need inla.tmarginal() here, as we actually look directly at the linear predictor

#Find mean of linear predictor from marginal

#Find Standard Deviation of linear predictor from marginal


#Calculate CRPS for each individual county at time t*: using mean and std of linear predictor
#i.e. use crpsNormal like crpsNormal(x = log(rate), mu = mean(eta), sig = sd(eta))


#Calculate absolute error for each individual county at time t*: using mean of linear predictor 















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

