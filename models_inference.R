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
load("proper_fitted.RData"); load("one_step_predictions.RData"); load("one_step_predictions_RW2.RData")



################################################################################
#Create latex table with results
#Format results data.frame for latex table

#Make one table containing all

cpo_waic_results.df <- data.frame(Model = rep(c("Improper 1 noInt",
                                              "Improper 1 typeI",
                                              "Improper 1 typeII",
                                              "Improper 1 typeIII",
                                              "Improper 1 typeIV",
                                              "Improper 2 noInt",
                                              "Improper 2 typeI",
                                              "Improper 2 typeII",
                                              "Improper 2 typeIII",
                                              "Improper 2 typeIV",
                                              "Proper noInt",
                                              "Proper onlyInt",
                                              "Proper full"), 3),
                                  model_choice = rep(c("1", "2", "3"), 13),
                                  value = 1:39)


cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 noInt" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 noInt" & cpo_waic_results.df$model_choice == "2"] = RW1_ICAR_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 noInt" & cpo_waic_results.df$model_choice == "3"] = time_RW1_ICAR

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeI" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_I_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeI" & cpo_waic_results.df$model_choice == "2"] = RW1_ICAR_I_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeI" & cpo_waic_results.df$model_choice == "3"] = time_RW1_ICAR_I

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeII" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_II_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeII" & cpo_waic_results.df$model_choice == "2"] = RW1_ICAR_II_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeII" & cpo_waic_results.df$model_choice == "3"] = time_RW1_ICAR_II

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeIII" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_III_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeIII" & cpo_waic_results.df$model_choice == "2"] = RW1_ICAR_III_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeIII" & cpo_waic_results.df$model_choice == "3"] = time_RW1_ICAR_III

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeIV" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW1_ICAR_IV_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeIV" & cpo_waic_results.df$model_choice == "2"] = RW1_ICAR_IV_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 1 typeIV" & cpo_waic_results.df$model_choice == "3"] = time_RW1_ICAR_IV

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 noInt" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 noInt" & cpo_waic_results.df$model_choice == "2"] = RW2_ICAR_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 noInt" & cpo_waic_results.df$model_choice == "3"] = time_RW2_ICAR

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeI" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_I_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeI" & cpo_waic_results.df$model_choice == "2"] = RW2_ICAR_I_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeI" & cpo_waic_results.df$model_choice == "3"] = time_RW2_ICAR_I

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeII" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_II_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeII" & cpo_waic_results.df$model_choice == "2"] = RW2_ICAR_II_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeII" & cpo_waic_results.df$model_choice == "3"] = time_RW2_ICAR_II

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeIII" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_III_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeIII" & cpo_waic_results.df$model_choice == "2"] = RW2_ICAR_III_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeIII" & cpo_waic_results.df$model_choice == "3"] = time_RW2_ICAR_III

cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeIV" & cpo_waic_results.df$model_choice == "1"] = mean(-log(RW2_ICAR_IV_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeIV" & cpo_waic_results.df$model_choice == "2"] = RW2_ICAR_IV_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Improper 2 typeIV" & cpo_waic_results.df$model_choice == "3"] = time_RW2_ICAR_IV

#Proper base results
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper noInt" & cpo_waic_results.df$model_choice == "1"] = mean(-log(proper_base_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper noInt" & cpo_waic_results.df$model_choice == "2"] = proper_base_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper noInt" & cpo_waic_results.df$model_choice == "3"] = time_proper_base

#Proper interaction results
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper onlyInt" & cpo_waic_results.df$model_choice == "1"] = mean(-log(proper_interaction_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper onlyInt" & cpo_waic_results.df$model_choice == "2"] = proper_interaction_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper onlyInt" & cpo_waic_results.df$model_choice == "3"] = time_proper_interaction

#Full proper results
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper full" & cpo_waic_results.df$model_choice == "1"] = mean(-log(proper_full_fit$cpo$cpo))
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper full" & cpo_waic_results.df$model_choice == "2"] = proper_full_fit$waic$waic
cpo_waic_results.df$value[cpo_waic_results.df$Model == "Proper full" & cpo_waic_results.df$model_choice == "3"] = time_proper_full * 60


#model = 1 -> base, model = 2 -> type I,..., model = 5 -> type IV
#RW = 1 -> RW1, RW = 2 -> RW2
#model_choice = 1 -> CPO, model_choice = 2 -> WAIC, model_choice = 3 -> Computational time (s)


#Make caption and label for latex table
caption = "HEIHEI"
label = "CPO-results"

#Make latex table
latex_tabular <- latexTable(tabular(
  Heading("Model")*RowFactor(Model, 
                        nopagebreak = "\\hline",
                        spacing = 0)~
    Heading()*Factor(model_choice, 
                     levelnames = c("CPO", "WAIC", "Computational time (s)"))*
    Heading()*value*Heading()*identity,
  data = cpo_waic_results.df),
  caption = caption,
  label = label
)
latex_tabular

#Save latex table
cat(latex_tabular, file = "cpo_table.tex")



################################################################################

############
##Plots of posteriors: save to pdf 7 by 3.2
plot_intercept(RW1_ICAR_fit, proper_base_fit)


#Plot the temporal random effect RW1 vs RW2: save to pdf 7.5 by 3.5
plot_temporal_effects_RW1_RW2(RW1_ICAR_fit, RW2_ICAR_fit, T)

#Plot temporal random effect ar1 + fixed (density fixed, ar1, ar1 + fixed)
#Save to pdf 10 by 3
plot_temporal_ar1(proper_base_fit)



#Plot spatial effects Proper vs improper 7.5 by 3.5
plot_spatial_effects(RW1_ICAR_fit,
                     proper_base_fit, 
                     ohio_map,
                     n)

#7.5 by 3.5
plot_spatial_std(RW1_ICAR_fit,
                 proper_base_fit, 
                 ohio_map,
                 n)



#####
#Plot the posterior hyperparameters (should it include RW2?)
#Plot improper temporal posterior hyperparameters: save to pdf 10 by 4
plot_improper_temporal_hyperparameters(RW1_ICAR_fit, RW2_ICAR_fit)


#Plot improper spatial posterior hyperparameters: 10 by 3.5
plot_improper_spatial_hyperparameters(RW1_ICAR_fit)


#Plot proper temporal hyperparameters: save to pdf 7.5 by 3.5
plot_proper_temporal_hyperparameter(proper_base_fit, proper_full_fit)


#Plot proper spatial hyperparameters: 7.5 by 3.5
plot_proper_spatial_hyperparameters(proper_base_fit, proper_full_fit)


########

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
#RW1, 7.5 by 3.5
plot_std_interactions_RW1(RW1_ICAR_I_fit,
                          RW1_ICAR_II_fit,
                          RW1_ICAR_III_fit,
                          RW1_ICAR_IV_fit)
#RW2
plot_std_interactions_RW2(RW2_ICAR_I_fit,
                          RW2_ICAR_II_fit,
                          RW2_ICAR_III_fit,
                          RW2_ICAR_IV_fit)


#Proper: save 11.25 by 3.5
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


counties = c(id_max_avg_rate, id_min_avg_rate, id_max_range_rate, 4)

#Produces plot of fitted vs true for
# first frame: county w. largest average rate
# second frame: county w. smallest average rate
# third frame: county w. largest range in rate
# fourth frame: a random county
# save pdf: 12 by 9 
select_county_timeseries(ohio_df,
                         RW1_ICAR_fit,
                         RW1_ICAR_IV_fit,
                         proper_interaction_fit,
                         counties,
                         n,
                         T)






#Do the same but for predicted values

predicted_vs_true_select_counties(base_predicted,
                                  IV_predicted,
                                  proper_interaction_predicted,
                                  pop_in_values_pred_on,
                                  values_predicted_on,
                                  counties)



#####
#Plot the heatmaps of some years (1968, 1973, 1978, 1983, 1988 maybe?) for some models
#Base model, RW1 type II, RW1 type IV, maybe proper ones


#Plot heatmaps of predicted values maybe???

years_to_plot = c(1979, 1982, 1985, 1988)
years_to_plot_id_pred = c(1, 4, 7, 10)

actual = ohio_df

#Extract true values for years_to_plot
if(actual$year[1]==1){
  actual$year = actual$year + (1968 - 1)
}

actual <- actual[actual$year %in% years_to_plot, ]

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

p1 <- case_count_plot_1_year(actual_n_map, years_to_plot[1], hardcoded_bins, title = "True rate: 1979")
p2 <- case_count_plot_1_year(actual_n_map, years_to_plot[2], hardcoded_bins, title = "True rate: 1982")
p3 <- case_count_plot_1_year(actual_n_map, years_to_plot[3], hardcoded_bins, title = "True rate: 1985")
p4 <- case_count_plot_1_year(actual_n_map, years_to_plot[4], hardcoded_bins, title = "True rate: 1988")

#Save as 11.25 by 3.5
ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1,
          common.legend = TRUE, legend = "right") 


#Calculate predicted values for Improper_1_typeIV predictions for these years



#Find predicted values for Improper_1_typeIV
predicted_IV_1979 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(IV_predicted[1, i])
  predicted_IV_1979[i] <- pred_rate * 1E5
}

predicted_IV_1982 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(IV_predicted[4, i])
  predicted_IV_1982[i] <- pred_rate * 1E5
}

predicted_IV_1985 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(IV_predicted[7, i])
  predicted_IV_1985[i] <- pred_rate * 1E5
}

predicted_IV_1988 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(IV_predicted[10, i])
  predicted_IV_1988[i] <- pred_rate * 1E5
}

actual_n_map$rate[actual_n_map$year == 1979] = predicted_IV_1979
actual_n_map$rate[actual_n_map$year == 1982] = predicted_IV_1982
actual_n_map$rate[actual_n_map$year == 1985] = predicted_IV_1985
actual_n_map$rate[actual_n_map$year == 1988] = predicted_IV_1988


plt1 <- case_count_plot_1_year(actual_n_map, years_to_plot[1], hardcoded_bins, title = "Improper_1_typeIV\n one-step prediction 1979")
plt2 <- case_count_plot_1_year(actual_n_map, years_to_plot[2], hardcoded_bins, title = "Improper_1_typeIV\n one-step prediction 1982")
plt3 <- case_count_plot_1_year(actual_n_map, years_to_plot[3], hardcoded_bins, title = "Improper_1_typeIV\n one-step prediction 1985")
plt4 <- case_count_plot_1_year(actual_n_map, years_to_plot[4], hardcoded_bins, title = "Improper_1_typeIV\n one-step prediction 1988")

ggarrange(plt1, plt2, plt3, plt4, ncol = 4, nrow = 1,
          common.legend = TRUE, legend = "right") 


#Find predicted values for Improper_1_typeIV
predicted_proper_1979 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(proper_interaction_predicted[1, i])
  predicted_proper_1979[i] <- pred_rate * 1E5
}

predicted_proper_1982 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(proper_interaction_predicted[4, i])
  predicted_proper_1982[i] <- pred_rate * 1E5
}

predicted_proper_1985 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(proper_interaction_predicted[7, i])
  predicted_proper_1985[i] <- pred_rate * 1E5
}

predicted_proper_1988 = rep(0, n)
for(i in 1:n){
  pred_rate = find_mean_marginal(proper_interaction_predicted[10, i])
  predicted_proper_1988[i] <- pred_rate * 1E5
}

actual_n_map$rate[actual_n_map$year == 1979] = predicted_proper_1979
actual_n_map$rate[actual_n_map$year == 1982] = predicted_proper_1982
actual_n_map$rate[actual_n_map$year == 1985] = predicted_proper_1985
actual_n_map$rate[actual_n_map$year == 1988] = predicted_proper_1988


plt1 <- case_count_plot_1_year(actual_n_map, years_to_plot[1], hardcoded_bins, title = "Proper_onlyInt\n one-step prediction 1979")
plt2 <- case_count_plot_1_year(actual_n_map, years_to_plot[2], hardcoded_bins, title = "Proper_onlyInt\n one-step prediction 1982")
plt3 <- case_count_plot_1_year(actual_n_map, years_to_plot[3], hardcoded_bins, title = "Proper_onlyInt\n one-step prediction 1985")
plt4 <- case_count_plot_1_year(actual_n_map, years_to_plot[4], hardcoded_bins, title = "Proper_onlyInt\n one-step prediction 1988")

ggarrange(plt1, plt2, plt3, plt4, ncol = 4, nrow = 1,
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


#####
#Fitted values of best model against true values 
plot_fitted_vs_actual_together(ohio_df, RW1_ICAR_fit,
                               RW1_ICAR_II_fit, n, T)

############
##One-step predictor ...

years_predicted_on <- 12:21

#Extract deaths from years which we predicted on
values_predicted_on <- ohio_df$deaths[ohio_df$year %in% years_predicted_on]

#Extract population in years predicted on
pop_in_values_pred_on <- ohio_df$pop_at_risk[ohio_df$year %in% years_predicted_on]




#Find MAE and RMSE
improper_1_noInt_MAE_RMSE = find_MAE_RMSE_all_years(base_predicted, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)


#Find IS
improper_1_IS <- find_IS_all(base_predicted, pop_in_values_pred_on, values_predicted_on)



improper_1_typeI_MAE_RMSE = find_MAE_RMSE_all_years(I_predicted, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)

improper_1_typeI_IS <- find_IS_all(I_predicted, pop_in_values_pred_on, values_predicted_on)


improper_1_typeII_MAE_RMSE = find_MAE_RMSE_all_years(II_predicted, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)


improper_1_typeII_IS <- find_IS_all(II_predicted, pop_in_values_pred_on, values_predicted_on)



improper_1_typeIII_MAE_RMSE = find_MAE_RMSE_all_years(III_predicted, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)

improper_1_typeIII_IS <- find_IS_all(III_predicted, pop_in_values_pred_on, values_predicted_on)

improper_1_typeIV_MAE_RMSE = find_MAE_RMSE_all_years(IV_predicted, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)

improper_1_typeIV_IS <- find_IS_all(IV_predicted, pop_in_values_pred_on, values_predicted_on)

#Find CRPS and AE for improper models with RW2
improper_2_noInt_MAE_RMSE = find_MAE_RMSE_all_years(improper_2_noInt_pred, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)

improper_2_IS <- find_IS_all(improper_2_noInt_pred, pop_in_values_pred_on, values_predicted_on)

improper_2_typeI_MAE_RMSE = find_MAE_RMSE_all_years(improper_2_I_pred, 
                                                  pop_in_values_pred_on,
                                                  values_predicted_on)

improper_2_typeI_IS <- find_IS_all(improper_2_I_pred, pop_in_values_pred_on, values_predicted_on)

improper_2_typeII_MAE_RMSE = find_MAE_RMSE_all_years(improper_2_II_pred, 
                                                   pop_in_values_pred_on,
                                                   values_predicted_on)

improper_2_typeII_IS <- find_IS_all(improper_2_II_pred, pop_in_values_pred_on, values_predicted_on)

improper_2_typeIII_MAE_RMSE = find_MAE_RMSE_all_years(improper_2_III_pred, 
                                                    pop_in_values_pred_on,
                                                    values_predicted_on)

improper_2_typeIII_IS <- find_IS_all(improper_2_III_pred, pop_in_values_pred_on, values_predicted_on)

improper_2_typeIV_MAE_RMSE = find_MAE_RMSE_all_years(improper_2_IV_pred, 
                                                   pop_in_values_pred_on,
                                                   values_predicted_on)

improper_2_typeIV_IS <- find_IS_all(improper_2_IV_pred, pop_in_values_pred_on, values_predicted_on)


#Find CRPS and AE for proper models
proper_noInt_MAE_RMSE <- find_MAE_RMSE_all_years(proper_base_predicted,
                                               pop_in_values_pred_on,
                                               values_predicted_on)

proper_noInt_IS <- find_IS_all(proper_base_predicted, pop_in_values_pred_on, values_predicted_on)

proper_interaction_MAE_RMSE <- find_MAE_RMSE_all_years(proper_interaction_predicted,
                                                     pop_in_values_pred_on,
                                                     values_predicted_on)

proper_interaction_IS <- find_IS_all(proper_interaction_predicted, pop_in_values_pred_on, values_predicted_on)

proper_full_MAE_RMSE <- find_MAE_RMSE_all_years(proper_full_predicted,
                                              pop_in_values_pred_on,
                                              values_predicted_on)

proper_full_IS <- find_IS_all(proper_full_predicted, pop_in_values_pred_on, values_predicted_on)



#Create a table with calculated AE and CRPS values

IS.df <- data.frame(Model = rep(c("Improper 1 noInt",
                                  "Improper 1 typeI",
                                  "Improper 1 typeII",
                                  "Improper 1 typeIII",
                                  "Improper 1 typeIV",
                                  "Improper 2 noInt",
                                  "Improper 2 typeI",
                                  "Improper 2 typeII",
                                  "Improper 2 typeIII",
                                  "Improper 2 typeIV",
                                  "proper noInt",
                                  "proper onlyInt",
                                  "proper full"), 3),
                      ae_crps = c(rep("MAE", 13), 
                                  rep("RMSE", 13), 
                                  rep("IS", 13)),
                      value = 1:39)


IS.df$value[IS.df$Model == "Improper 1 noInt" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_1_noInt_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 1 noInt" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_1_noInt_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 1 noInt" &
              IS.df$ae_crps == "IS"] = round(improper_1_IS, 2)


IS.df$value[IS.df$Model == "Improper 1 typeI" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_1_typeI_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 1 typeI" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_1_typeI_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 1 typeI" &
              IS.df$ae_crps == "IS"] = round(improper_1_typeI_IS, 2)



IS.df$value[IS.df$Model == "Improper 1 typeII" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_1_typeII_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 1 typeII" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_1_typeII_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 1 typeII" &
              IS.df$ae_crps == "IS"] = round(improper_1_typeII_IS, 2)


IS.df$value[IS.df$Model == "Improper 1 typeIII" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_1_typeIII_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 1 typeIII" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_1_typeIII_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 1 typeIII" &
              IS.df$ae_crps == "IS"] = round(improper_1_typeIII_IS, 2)



IS.df$value[IS.df$Model == "Improper 1 typeIV" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_1_typeIV_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 1 typeIV" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_1_typeIV_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 1 typeIV" &
              IS.df$ae_crps == "IS"] = round(improper_1_typeIV_IS, 2)




IS.df$value[IS.df$Model == "Improper 2 noInt" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_2_noInt_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 2 noInt" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_2_noInt_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 2 noInt" &
              IS.df$ae_crps == "IS"] = round(improper_2_IS, 2)


IS.df$value[IS.df$Model == "Improper 2 typeI" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_2_typeI_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 2 typeI" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_2_typeI_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 2 typeI" &
              IS.df$ae_crps == "IS"] = round(improper_2_typeI_IS, 2)



IS.df$value[IS.df$Model == "Improper 2 typeII" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_2_typeII_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 2 typeII" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_2_typeII_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 2 typeII" &
              IS.df$ae_crps == "IS"] = round(improper_2_typeII_IS, 2)


IS.df$value[IS.df$Model == "Improper 2 typeIII" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_2_typeIII_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 2 typeIII" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_2_typeIII_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 2 typeIII" &
              IS.df$ae_crps == "IS"] = round(improper_2_typeIII_IS, 2)



IS.df$value[IS.df$Model == "Improper 2 typeIV" &
              IS.df$ae_crps == "MAE"] = round(mean(improper_2_typeIV_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "Improper 2 typeIV" &
              IS.df$ae_crps == "RMSE"] = round(mean(improper_2_typeIV_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "Improper 2 typeIV" &
              IS.df$ae_crps == "IS"] = round(improper_2_typeIV_IS, 2)





IS.df$value[IS.df$Model == "proper noInt" &
              IS.df$ae_crps == "MAE"] = round(mean(proper_noInt_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "proper noInt" &
              IS.df$ae_crps == "RMSE"] = round(mean(proper_noInt_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "proper noInt" &
              IS.df$ae_crps == "IS"] = round(proper_noInt_IS, 2)


IS.df$value[IS.df$Model == "proper onlyInt" &
              IS.df$ae_crps == "MAE"] = round(mean(proper_interaction_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "proper onlyInt" &
              IS.df$ae_crps == "RMSE"] = round(mean(proper_interaction_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "proper onlyInt" &
              IS.df$ae_crps == "IS"] = round(proper_interaction_IS, 2)


IS.df$value[IS.df$Model == "proper full" &
              IS.df$ae_crps == "MAE"] = round(mean(proper_full_MAE_RMSE$MAE), 2)
IS.df$value[IS.df$Model == "proper full" &
              IS.df$ae_crps == "RMSE"] = round(mean(proper_full_MAE_RMSE$RMSE), 2)

IS.df$value[IS.df$Model == "proper full" &
              IS.df$ae_crps == "IS"] = round(proper_full_IS, 2)



#Make caption and label for latex table
caption_IS = "HEIHEI"
label_IS = "mae-rmse-is"


#Make table for proper results

#Make latex table
latex_tabular_IS <- latexTable(tabular(
  Heading("Model")*RowFactor(Model,
                             nopagebreak = "\\hline",
                             spacing = 0)~
    Heading()*Factor(ae_crps)*
    Heading()*value*Heading()*identity,
  data = IS.df),
  caption = caption_IS,
  label = label_IS
)
latex_tabular_IS

#Save latex table
cat(latex_tabular_IS, file = "IS.tex")











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

