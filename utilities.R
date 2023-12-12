#Constraint maker function
constraints_maker <- function(type = NULL, n = NULL, t = NULL,
                              rw = "RW1", prec_matrix = NULL){
  #Type: specifies what interaction type and hence what type of constraint is desired
  #n: specifies number of areas
  #t: specifies number of time points
  #rw: specifies if temporal random effects follows a RW1 or RW2
  #prec_matrix: Precision matrix, used to define constraints using eigenvectors (only for RW2)
  
  if(rw == "RW1"){ #Define constraints for RW(1)
    if(type == "II"){
      #For a type II interaction, there is a RW(1) over each interaction (assuming RW1)
      #Hence for each county, constrain the RW to sum-to-zero
      A <- matrix(0, nrow = n, ncol = n * t)
      for (i in 1:(n - 1)) {
        A[i, which((1:(n * t))%%n == i)] <- 1
      }
      A[n, which((1:(n * t))%%n == 0)] <- 1
      
    } else if(type == "III"){
      #For a type III interaction, there is a indep. ICAR at each time point
      #Need the ICAR at each time point to sum-to-zero
      A <- matrix(0, nrow = t, ncol = n * t)
      for (i in 1:t) {
        # The ICAR at each time point needs to sum to 0
        A[i, ((i - 1) * n + 1):(i * n)] <- 1
      }
      
    } else if(type == "IV"){
      #For a type IV interaction, we have to do both sum-to-zero
      #over each RW on each county, and for each time point sum-to-zero
      #over each ICAR
      time_constr <- matrix(0, nrow = n, ncol = n * t)
      for (i in 1:(n - 1)) {
        time_constr[i, which((1:(n * t))%%n == i)] <- 1
      }
      time_constr[n, which((1:(n * t))%%n == 0)] <- 1
      
      space_constr <- matrix(0, nrow = t-1, ncol = n * t)
      for (i in 1:(t-1)) { 
        space_constr[i, ((i - 1) * n + 1):(i * n)] <- 1
      }
      
      A <- rbind(time_constr, space_constr)
    }
  } else{ #Define constraints with RW2
    if(type == "II"){
      #For a type II interaction, there is a RW(2) over each interaction (assuming RW2)
      eigens <- eigen(prec_matrix)
      
      #Extract 2n last eigenvectors corresponding to eigenvalues=0
      A <- t(eigens$vectors[ ,(nrow(eigens$vectors) - 2 * n + 1):nrow(eigens$vectors)])
      
    } else if(type == "III"){
      #For a type III interaction, there is a indep. ICAR at each time point
      #Need the ICAR at each time point to sum-to-zero
      A <- matrix(0, nrow = t, ncol = n * t)
      for (i in 1:t) {
        # The ICAR at each time point needs to sum to 0
        A[i, ((i - 1) * n + 1):(i * n)] <- 1
      }
      
    } else if(type == "IV"){
      #Calculate eigenvectors
      eigens <- eigen(prec_matrix)
      
      #Extract 2n + T - 2 last eigenvectors corresponding to eigenvalues=0
      A <- t(eigens$vectors[ ,(nrow(eigens$vectors) - 2 * n - t + 3):nrow(eigens$vectors)])
    }
  }
  
  #Get constraints in INLA format
  constr.st <- list(A = A, e = rep(0, dim(A)[1]))
  return(constr.st)
}

#Function printing some key numbers and stuff
print_cpo_etc <- function(fitted_model, time_obj){
  print(c("Number of constraints: ",toString(fitted_model$misc$configs$constr$nc)))
  print(c("- sum(log(CPO)):", toString(round(-sum(log(fitted_model$cpo$cpo)), digits = 0))))
  print(c("WAIC: ", toString(round(fitted_model$waic$waic, digits = 0))))
  print(c("CPO failiure, should have 1 unique value",
          toString(length(unique(fitted_model$cpo$failure))))) #-> soinla.cpo(res)
  print(c("And that value is (should be 0): ", toString(fitted_model$cpo$failure[1])))
  print(c("CPU: ", toString(summary(fitted_model)$cpu.used)))
  print(c("Time to compute", toString(round(time_obj, 3))))
}

#Function to sort proper models so that data appears in same order as for improper
sort_proper_fitted <- function(proper_fitted, n, time){
  sorted_proper_fitted <- proper_fitted
  for(t in 1:time){
    sorted_proper_fitted[((t-1)*n + 1):(t*n), ] =  proper_fitted[seq(t, n*time, by = time), ]
  }
  return(sorted_proper_fitted)
}

sort_proper_fitted_2 <- function(proper_fitted, n, time){
  sorted_proper_fitted <- proper_fitted
  for(t in 1:time){
    sorted_proper_fitted[((t-1)*n + 1):(t*n)] =  proper_fitted[seq(t, n*time, by = time)]
  }
  return(sorted_proper_fitted)
}

#Function to find mean fitted linear predictor marginal for singular instance
find_mean_marginal <- function(marginal){
  return(mean(marginal[[1]][, 1]))
}

#Function to find standard deviation of fitted linear predictor marginal for singular instance
find_var_marginal <- function(marginal){
  return(var(marginal[[1]][, 1]))
}

#Function used to go from linear predictor marginal to fitted rate marginal
my_inla_t_marginal <- function(prediction_marginal){
  #a function to use inla.tmarginal on several values at once
  return(inla.tmarginal(function(x){exp(x)}, prediction_marginal))
}

square_error <- function(x, mu = 0){
  return((x - mu)^2)
}

abs_error <- function(x, mu = 0){
  return(abs(x - mu))
}

MAE_and_MSE_one_year <- function(marginals, pop, true_values, year){
  
  abs_error = rep(0, n); square_error = rep(0, n)
  
  mean_marg = rep(0, n); mean_predicted = rep(0, n)

  
  for(i in 1:n){
    mean_marg[i] = find_mean_marginal(marginals[year, i])
    mean_predicted[i] = pop[((year - 1) * n + i)] * mean_marg[i]
    
    abs_error[i] = abs_error(true_values[((year - 1) * n + i)],
                             mu = mean_predicted[i])
    
    square_error[i] = square_error(true_values[((year - 1) * n + i)],
                                   mu = mean_predicted[i])
    
    
  }
  return(list(ae = mean(abs_error),
              mse = mean(square_error)))
}


find_MAE_RMSE_all_years <- function(marginals, pop, true){
  AE = rep(0, 10)
  MSE = rep(0, 10)
  for(t in 1:10){
    temp   = MAE_and_MSE_one_year(marginals, pop, true, t)
    AE[t]  = mean(temp$ae)
    MSE[t] = mean(temp$mse)
  }
  MAE = mean(AE)
  RMSE = sqrt(mean(MSE))
  return(list(MAE = MAE, 
              RMSE = RMSE))
}


#Sample from posterior marginal distribution of lambda using inla.rmarginal,
#Then sample mortality count O_it* from poisson w. distribution n_{it}lambda sampled ones for each sample of lambda
#Use the sampled counts to produce upper and lower quantiles
find_u_l_single_pred <- function(marginal_lambda, population, n_samples = 5000){
  
  #Sample lambda from marginal
  samples <- inla.rmarginal(n_samples, marginal_lambda)
  
  #scale by population
  samples_offset <- population * samples
  
  #Sample from poisson using samples_offset
  poisson_samples = rep(0, n_samples)
  for(i in 1:n_samples){
    poisson_samples[i] <- rpois(1, samples_offset[i])
  }
  
  #Calculate upper and lower quantile
  u = as.numeric(quantile(poisson_samples, 0.975))
  l = as.numeric(quantile(poisson_samples, 0.025))
  
  return(list(l = l, u = u))
}

#test = find_u_l_single_pred(base_predicted[1, 1][[1]], pop_in_values_pred_on[1])

#Interval-scores
find_IS_one_obs <- function(l, u, true_value){
  return((u - l) + 2/0.05 * (l - true_value) * (true_value < l) + 2/0.05 * (true_value - u) * (true_value > u))
}

#lu = find_u_l_single_pred(base_predicted[1, 1][[1]], pop_in_values_pred_on[1])
#test = find_IS_one_obs(lu$l, lu$u, values_predicted_on[1])

find_IS_one_year <- function(marginals, population, true_values, year){
  IS = rep(0, n)
  
  for(i in 1:n){
    
    temp_lu = find_u_l_single_pred(marginals[year, i][[1]], population[((year - 1) * n + i)])
    IS[i] = find_IS_one_obs(temp_lu$l, temp_lu$u, true_values[((year - 1) * n + i)])
    
  }
  return(mean(IS))
}

#test = find_IS_one_year(base_predicted, pop_in_values_pred_on, values_predicted_on, 1)

find_IS_all <- function(marginals, population, true_values){
  #Iterate over years
  IS = rep(0, 10)
  for(t in 1:10){
    temp   = find_IS_one_year(marginals, population, true_values, t)
    IS[t]  = temp
  }
  return(mean(IS))
}

#test = find_IS_all(base_predicted, pop_in_values_pred_on, values_predicted_on)




#Plot the intercept
plot_intercept <- function(improper_base, proper_base){
  #Format intercept for ggplot
  improper.df <- data.frame(x_axis = improper_base$marginals.fixed$`(Intercept)`[, 1],
                            y_axis = improper_base$marginals.fixed$`(Intercept)`[, 2])
  
  proper.df <- data.frame(x_axis = proper_base$marginals.fixed$`(Intercept)`[, 1],
                          y_axis = proper_base$marginals.fixed$`(Intercept)`[, 2])
  
  improper_plot <- ggplot(data = improper.df, 
                          aes(x = x_axis, y = y_axis)) +
                  geom_area(fill = "#F8766D", alpha = 0.8) + #
                  theme_bw() +
                  theme(axis.title=element_text(size=14,face="bold")) +
                  xlab(expression(mu)) + ylab(expression(f(mu))) + 
                  ggtitle("Intercept improper models")
                          
  
  
  proper_plot <- ggplot(data = proper.df, 
                        aes(x = x_axis, y = y_axis)) +
                geom_area(fill = "#00BFC4", alpha = 0.8) + #
                theme_bw() + 
                theme(axis.title=element_text(size=14,face="bold")) +
                xlab(expression(mu)) + ylab(expression(f(mu))) + 
                ggtitle("Intercept proper models")
  
  
  ggarrange(improper_plot, proper_plot,
            ncol = 2, nrow = 1)
  
}

#Plot temporal effect RW1 vs RW2
plot_temporal_effects_RW1_RW2 <- function(fitted_RW1, fitted_RW2, T){
  #Format temporal random effects for ggplot
  years <- 1968:1988
  temporal_RW1.df <- data.frame(years = years)
  temporal_RW1.df$lower_quant <- fitted_RW1$summary.random$year[(T + 1):(2 * T), 4]
  temporal_RW1.df$median <- fitted_RW1$summary.random$year[(T + 1):(2 * T), 5]
  temporal_RW1.df$upper_quant <- fitted_RW1$summary.random$year[(T + 1):(2 * T), 6]
  
  temporal_RW2.df <- data.frame(years = years)
  temporal_RW2.df$lower_quant <- fitted_RW2$summary.random$year[(T + 1):(2 * T), 4]
  temporal_RW2.df$median <- fitted_RW2$summary.random$year[(T + 1):(2 * T), 5]
  temporal_RW2.df$upper_quant <- fitted_RW2$summary.random$year[(T + 1):(2 * T), 6]
  
  
  temporal_RW1 <- ggplot(data = temporal_RW1.df, aes(years, median)) + 
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line() + 
    geom_line(data = temporal_RW1.df, aes(years, lower_quant), linetype = "dashed") + 
    geom_line(data = temporal_RW1.df, aes(years, upper_quant), linetype = "dashed") + 
    xlab("year") + ylab(expression(alpha[t])) +
    ggtitle("RW1")
  
  temporal_RW2 <- ggplot(data = temporal_RW2.df, aes(years, median)) + 
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    geom_line() + 
    geom_line(data = temporal_RW2.df, aes(years, lower_quant), linetype = "dashed") + 
    geom_line(data = temporal_RW2.df, aes(years, upper_quant), linetype = "dashed") + 
    xlab("year") + ylab(expression(alpha[t])) + 
    ggtitle("RW2")
  
  ggarrange(temporal_RW1, temporal_RW2,
            ncol = 2, nrow = 1)
}


#Plot temporal effect ar1 + fixed
plot_temporal_ar1 <- function(proper_base){
  
  #proper_base_fit$marginals.fixed$year
  fixed.df <- data.frame(x_axis = proper_base$marginals.fixed$year[, 1],
                         y_axis = proper_base$marginals.fixed$year[, 2])
  
  
  fixed_plot <- ggplot(data = fixed.df, 
                       aes(x = x_axis, y = y_axis)) +
               geom_area(fill = "#999999", alpha = 0.8) +
               theme_bw() +
               theme(axis.title=element_text(size=14)) +
               xlab(expression(beta)) + ylab(expression(f(beta))) + 
               ggtitle("Density fixed effect")
  
  years <- 1968:1988
  
  ar1.df <- data.frame(year = years,
                       lower_quant = proper_base$summary.random$year.copy[, 4],
                       median = proper_base$summary.random$year.copy[, 5],
                       upper_quant = proper_base$summary.random$year.copy[, 6])
  
  
  ar1_plot <- ggplot(data = ar1.df, aes(year, median)) +
              theme_bw() + 
              theme(axis.title=element_text(size=14)) +
              geom_line() + 
              geom_line(data = ar1.df, aes(years, lower_quant), linetype = "dashed") + 
              geom_line(data = ar1.df, aes(years, upper_quant), linetype = "dashed") + 
              xlab("year") + ylab(expression(alpha[t])) +
              ggtitle("AR1")
  
  years_standardized = 1:21
  beta_mean = proper_base$summary.fixed$mean[2]
  ar1_fixed_median = ar1.df$median + beta_mean * years_standardized
  lower_quant = ar1.df$lower_quant + beta_mean * years_standardized
  upper_quant = ar1.df$upper_quant + beta_mean * years_standardized
  
  ar1_fixed.df = data.frame(year = years,
                            median = ar1_fixed_median,
                            lower_quant = lower_quant,
                            upper_quant = upper_quant)
  
  ar1_fixed_plot <- ggplot(data = ar1_fixed.df, aes(year, median)) + 
                    theme_bw() + 
                    theme(axis.title=element_text(size=14)) +
                    geom_line() + 
                    geom_line(data = ar1_fixed.df,
                              aes(years, lower_quant), linetype = "dashed") + 
                    geom_line(data = ar1_fixed.df,
                              aes(years, upper_quant), linetype = "dashed") + 
                    xlab("year") + ylab(expression(alpha[t]+beta*t)) + 
                    ggtitle("AR1 and fixed effect")
    
    
  #ggarrange(fixed_plot, ar1_plot, ar1_fixed_plot, 
  #          ncol = 3, nrow = 1)
  
  
  
  ggarrange(fixed_plot, ar1_plot,                                                 # First row with scatter plot
            ar1_fixed_plot, 
            nrow = 2,  ncol = 2) 
}



#Plot posterior intercept, temporal effects and spatial effects
plot_spatial_effects <- function(improper,
                                 proper,
                                 map,
                                 n){
  #Function that produces four plots: The posterior intercept, 
  #posterior structured temporal effect along with 2.5% and 97.5% quantiles
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  
  
  #scale_improper = scale_col[c(3, 10, 13, 18, 21, 24, 27, 30)]
  #scale_proper = scale_col[c(3,10,13,18,21,24,27,30)]
  
  
  scale_improper = scale_col[seq(1, 30, length.out = 10)]
  scale_proper = scale_col[seq(1, 30, length.out = 10)]
  
  improper_mean <- improper$summary.random$county$mean[(n+1):(2*n)]
  proper_mean <- proper$summary.random$county$mean
  
  temp_ohio_map <- map[ ,c("geometry", "NAME")]
  temp_ohio_map$improper_mean <- improper_mean
  temp_ohio_map$proper_mean <- proper_mean
  
  #Hardcoded bins, so as both heatmaps have same scale
  min_val_improper <- min(improper_mean)
  max_val_imporper <- max(improper_mean)
  min_val_proper <- min(proper_mean)
  max_val_proper <- max(proper_mean)
  
  hardcoded_bins_improper = round(seq(min_val_improper,
                                      max_val_imporper, length.out = 8), 2)
  
  
  #hardcoded_bins_improper = round(seq(0, 1, length.out = 10), 2)
  
  hardcoded_bins_proper = round(seq(min_val_proper,
                                    max_val_proper, length.out = 8), 2)
  
  
  #Set the theme to minimal
  theme_set(theme(panel.background = element_blank()))
  improper_mean_plot <- ggplot(data = temp_ohio_map) + 
                        geom_sf(aes(fill = improper_mean), 
                                alpha = 1,
                                color="black") + ggtitle("Mean ICAR") +
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
                          palette = function(x) c(scale_improper),
                          labels = function(x){x},
                          breaks = hardcoded_bins_improper,
                          limits = c(-3, 3),
                          guide = "colorscale")
  
  
  proper_mean_plot <- ggplot(data = temp_ohio_map) + 
                      geom_sf(aes(fill = proper_mean), 
                              alpha = 1,
                              color="black") + ggtitle("Mean CAR") +
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
                        palette = function(x) c(scale_proper),
                        labels = function(x){x},
                        breaks = hardcoded_bins_proper,
                        limits = c(-3, 3),
                        guide = "colorscale")
  
  ggarrange(improper_mean_plot, proper_mean_plot,
            ncol = 2, nrow = 1,
            common.legend = FALSE)
}

plot_spatial_std <- function(improper,
                             proper,
                             map,
                             n){
  
  #Function that produces four plots: The posterior intercept, 
  #posterior structured temporal effect along with 2.5% and 97.5% quantiles
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  
  #scale_1 = scale_col[c(3,10,13,17,21,24,27,30)] #Select color scale to be more red
  scale_1 = scale_col[seq(3, 30, length.out = 8)]
  
  
  improper_sd <- improper$summary.random$county$sd[(n+1):(2*n)]
  proper_sd <- proper$summary.random$county$sd
  
  temp_ohio_map <- map[ ,c("geometry", "NAME")]
  temp_ohio_map$improper_sd <- improper_sd
  temp_ohio_map$proper_sd <- proper_sd
  
  
  #Hardcoded bins
  #Hardcoded bins, so as both heatmaps have same scale
  min_val_improper <- min(improper_sd)
  max_val_imporper <- max(improper_sd)
  min_val_proper <- min(proper_sd)
  max_val_proper <- max(proper_sd)
  
  hardcoded_bins_improper = round(seq(min_val_improper - 0.01,
                                      max_val_imporper + 0.01, length.out = 8), 2)
  
  
  #hardcoded_bins_improper = round(seq(0, 1, length.out = 10), 2)
  
  hardcoded_bins_proper = round(seq(min_val_proper,
                                    max_val_proper, length.out = 8), 2)
  
  
  #Set the theme to minimal
  theme_set(theme(panel.background = element_blank()))
  improper_sd_plot <- ggplot(data = temp_ohio_map) + 
    geom_sf(aes(fill = improper_sd), 
            alpha = 1,
            color="black") + ggtitle("Standard deviation ICAR") +
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
      breaks = hardcoded_bins_improper,
      limits = c(0, 3),
      guide = "colorscale")
  
  proper_sd_plot <- ggplot(data = temp_ohio_map) + 
    geom_sf(aes(fill = proper_sd), 
            alpha = 1,
            color="black") + ggtitle("Standard deviation CAR") +
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
      breaks = hardcoded_bins_proper,
      limits = c(0, 3),
      guide = "colorscale")

  
  ggarrange(improper_sd_plot, proper_sd_plot,
            ncol = 2, nrow = 1,
            common.legend = FALSE)
  
}









#Sample posterior to get phi times variance of county???
hyperpar_sampler <- function(fitted_model,
                             n_samples = 10000,
                             temporal = TRUE){
  
  #Matrix where each row is a sample, and columns is which hyperparameter it is a sample of
  posterior_samples = inla.hyperpar.sample(n_samples, fitted_model) #Sample the posterior
  
  if(temporal){ #We want the temporal hyperparameters
    #Precision for year; Phi for year, col 1 and 2
    #Want to find standard deviation explained by alpha_t
    std_explained = rep(0, n_samples)
    for(i in 1:n_samples){
      std_explained[i] = sqrt(posterior_samples[i, 2] * as.numeric(1/posterior_samples[i, 1]))
    }
    
    #Get variance explained by structured effect, i.e. product of marginal variance and Phi for year
    return(std_explained)
    
    
  } else{ #We want the spatial hyperparameters
    #Precision for county, Phi for county col 3 and 4 respectively
    std_explained = rep(0, n_samples)
    
    for(i in 1:n_samples){
      std_explained[i] = sqrt(posterior_samples[i, 4] * as.numeric(1/posterior_samples[i, 3]))
    }
    
    #Get variance explained by structured effect, i.e. product of marginal variance and Phi for year
    return(std_explained)
  }
  
}




#Plot posterior hyperparameters of temporal random effects
plot_improper_temporal_hyperparameters <- function(fitted_RW1, fitted_RW2){
  #inla.tmarginal to transform from precision to standard deviation
  std_year_RW1 <- inla.tmarginal(function(x) sqrt(1/x),
                                 fitted_RW1$marginals.hyperpar$`Precision for year`)
  
  std_year_RW2 <- inla.tmarginal(function(x) sqrt(1/x),
                                 fitted_RW2$marginals.hyperpar$`Precision for year`)
  
  
  #Format for ggplot
  std_temporal_df <- data.frame(x_axis = c(std_year_RW1[, 1], std_year_RW2[, 1]),
                                y_axis = c(std_year_RW1[, 2], std_year_RW2[, 2]),
                                type = c(rep("RW1", length(std_year_RW1[, 1])),
                                         rep("RW2", length(std_year_RW2[, 1]))))
  
  
  phi_temporal_df <- data.frame(x_axis = c(fitted_RW1$marginals.hyperpar$`Phi for year`[, 1],
                                           fitted_RW2$marginals.hyperpar$`Phi for year`[, 1]),
                                y_axis = c(fitted_RW1$marginals.hyperpar$`Phi for year`[, 2],
                                           fitted_RW2$marginals.hyperpar$`Phi for year`[, 2]),
                                type = c(rep("RW1", 
                                             length(fitted_RW1$marginals.hyperpar$`Phi for year`[, 1])),
                                         rep("RW2", 
                                             length(fitted_RW2$marginals.hyperpar$`Phi for year`[, 1]))))
  
  std_temporal_plot <- ggplot(data=std_temporal_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(TeX("$f(*)$")) #expression(f(sigma)) 
  
  phi_temporal_plot <- ggplot(data=phi_temporal_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=6.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(phi)) + ylab("") #expression(f(lambda))
  
  
  
  #Extract how much of the standard deviation is captured
  samples_rw1 = hyperpar_sampler(fitted_RW1)
  samples_rw2 = hyperpar_sampler(fitted_RW2)
  
  sigma_sqrtPhi.df <- data.frame(x_axis = c(samples_rw1,
                                            samples_rw2),
                                 type = c(rep("RW1", length(samples_rw1)),
                                          rep("RW2", length(samples_rw2))))

  sigma_phi.plot <- ggplot(data=sigma_sqrtPhi.df) + 
                            stat_density(aes(x=x_axis, group=type, fill=type),
                                         adjust=1.5, alpha=.8, position = "identity") +
                            theme_bw() +
                            theme(axis.title=element_text(size=14)) +
                            labs(fill = NULL) +
                            xlab(TeX("$\\sigma\\sqrt{\\phi}$")) + #"$\\left(\\frac{\\phi}{\\tau}\\right)^{1/2}$"
                            ylab("")
  
  plt <- ggarrange(std_temporal_plot, phi_temporal_plot, sigma_phi.plot,
                    ncol = 3, nrow = 1,
                    common.legend = T, legend = "top")
  
  annotate_figure(plt, top = text_grob("Temporal hyperparameters of improper models", 
                                       color = "black", size = 14))
  
  
}


#Plot posterior hyperparameters of spatial random effects
plot_improper_spatial_hyperparameters <- function(fitted_RW1){
  std_county_RW1 <- inla.tmarginal(function(x) sqrt(1/x),
                                   fitted_RW1$marginals.hyperpar$`Precision for county`)
  
  std_spatial_df <- data.frame(x_axis = std_county_RW1[, 1],
                               y_axis = std_county_RW1[, 2])
  
  phi_spatial_df <- data.frame(x_axis = fitted_RW1$marginals.hyperpar$`Phi for county`[, 1],
                               y_axis = fitted_RW1$marginals.hyperpar$`Phi for county`[, 2])
  
  std_spatial_plot <- ggplot(data=std_spatial_df,
                             aes(x=x_axis)) +
    stat_density(fill = "#F8766D", adjust=1.5, alpha=.8,
                 position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(TeX("$f(*)$")) #expression(f(sigma))
  
  phi_spatial_plot <- ggplot(data=phi_spatial_df,
                             aes(x=x_axis)) +
    stat_density(fill = "#00BFC4", adjust=3, alpha=.8,
                 position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(phi)) + ylab("") #expression(f(lambda))
  
  
  #add w. 
  #Extract how much of the standard deviation is captured
  posterior_samples = hyperpar_sampler(fitted_RW1, temporal = FALSE)
  
  sigma_sqrtPhi.df <- data.frame(x_axis = posterior_samples)
  
  sigma_phi.plot <- ggplot(data=sigma_sqrtPhi.df,
                           aes(x = x_axis)) + 
                  stat_density(fill = "lightgreen",
                               adjust=1.5, alpha=.8, position = "identity") +
                  theme_bw() +
                  theme(axis.title=element_text(size=14)) +
                  labs(fill = NULL) +
                  xlab(TeX("$\\sigma\\sqrt{\\phi}$")) + 
                  ylab("")
  
  
  plt <- ggarrange(std_spatial_plot, phi_spatial_plot, sigma_phi.plot,
                   ncol = 3, nrow = 1)
  
  annotate_figure(plt, top = text_grob("Spatial hyperparameters of improper models", 
                  color = "black", size = 14))
}

#Plot posterior hyperparameters of proper temporal random effect
plot_proper_temporal_hyperparameter <- function(proper_base, proper_full){
  
  std_ar1_base <- inla.tmarginal(function(x) sqrt(1/x),
                                 proper_base$marginals.hyperpar$`Precision for year.copy`)
  
  std_ar1_full <- inla.tmarginal(function(x) sqrt(1/x),
                                 proper_full$marginals.hyperpar$`Precision for year.copy`)
  
  #Format for ggplot
  std_temporal_df <- data.frame(x_axis = c(std_ar1_base[, 1], std_ar1_full[, 1]),
                                y_axis = c(std_ar1_base[, 2], std_ar1_full[, 2]),
                                type = c(rep("Proper_noInt", length(std_ar1_base[, 1])),
                                         rep("Proper_full", length(std_ar1_full[, 1]))))
  
  
  rho_temporal_df <- data.frame(x_axis = c(proper_base$marginals.hyperpar$`Rho for year.copy`[, 1],
                                           proper_full$marginals.hyperpar$`Rho for year.copy`[, 1]),
                                y_axis = c(proper_base$marginals.hyperpar$`Rho for year.copy`[, 2],
                                           proper_full$marginals.hyperpar$`Rho for year.copy`[, 2]),
                                type = c(rep("Base", 
                                             length(proper_base$marginals.hyperpar$`Rho for year.copy`[, 1])),
                                         rep("Full", 
                                             length(proper_full$marginals.hyperpar$`Rho for year.copy`[, 1]))))
  
  std_temporal_plot <- ggplot(data=std_temporal_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=1, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) 
  
  rho_temporal_plot <- ggplot(data=rho_temporal_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=3, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho)) + ylab(expression(f(rho)))
  
  plt <- ggarrange(std_temporal_plot, rho_temporal_plot,
                  common.legend = T, legend = "top")
  
  annotate_figure(plt, top = text_grob("Temporal hyperparameters of proper models", 
                                       color = "black", size = 12))
  
}


#Plot posterior hyperparameters of proper spatial random effect
plot_proper_spatial_hyperparameters <- function(proper_base, proper_full){
  
  std_CAR_base <- inla.tmarginal(function(x) sqrt(1/x),
                                 proper_base$marginals.hyperpar$`Precision for county`)
  
  std_CAR_full <- inla.tmarginal(function(x) sqrt(1/x),
                                 proper_full$marginals.hyperpar$`Precision for county`)
  
  
  
  
  #Format for ggplot
  std_spatial_df <- data.frame(x_axis = c(std_CAR_base[, 1], std_CAR_full[, 1]),
                               y_axis = c(std_CAR_base[, 2], std_CAR_full[, 2]),
                               type = c(rep("Proper_noInt", length(std_CAR_base[, 1])),
                                        rep("Proper_full", length(std_CAR_full[, 1]))))
  
  
  lambda_spatial_df <- data.frame(x_axis = c(proper_base$marginals.hyperpar$`Lambda for county`[, 1],
                                             proper_full$marginals.hyperpar$`Lambda for county`[, 1]),
                                  y_axis = c(proper_base$marginals.hyperpar$`Lambda for county`[, 2],
                                             proper_full$marginals.hyperpar$`Lambda for county`[, 2]),
                                  type = c(rep("Base", 
                                               length(proper_base$marginals.hyperpar$`Lambda for county`[, 1])),
                                           rep("Full", 
                                               length(proper_full$marginals.hyperpar$`Lambda for county`[, 1]))))
  
  std_spatial_plot <- ggplot(data=std_spatial_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) 
  
  lambda_spatial_plot <- ggplot(data=lambda_spatial_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=2, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(lambda)) + ylab(expression(f(lambda)))
  
  plt <- ggarrange(std_spatial_plot, lambda_spatial_plot,
            common.legend = T, legend = "top")
  
  annotate_figure(plt, top = text_grob("Spatial hyperparameters of proper models", 
                                       color = "black", size = 12))
}




#Plot the interactions
plot_improper_interaction <- function(rw1_I, rw1_II, rw1_III, rw1_IV,
                             rw2_I, rw2_II, rw2_III, rw2_IV){
  x_axis = 1:length(rw1_I$summary.random$space.time[ , 4])
  #Format for ggplot
  rw1_I.df <- data.frame(x_axis = x_axis,
                         lower_quant = rw1_I$summary.random$space.time[ , 4], 
                         median = rw1_I$summary.random$space.time[ , 5],
                         upper_quant = rw1_I$summary.random$space.time[ , 6])
  
  rw1_II.df <- data.frame(x_axis = x_axis,
                         lower_quant = rw1_II$summary.random$space.time[ , 4], 
                         median = rw1_II$summary.random$space.time[ , 5],
                         upper_quant = rw1_II$summary.random$space.time[ , 6])
  
  rw1_III.df <- data.frame(x_axis = x_axis,
                         lower_quant = rw1_III$summary.random$space.time[ , 4], 
                         median = rw1_III$summary.random$space.time[ , 5],
                         upper_quant = rw1_III$summary.random$space.time[ , 6])
  
  rw1_IV.df <- data.frame(x_axis = x_axis,
                         lower_quant = rw1_IV$summary.random$space.time[ , 4], 
                         median = rw1_IV$summary.random$space.time[ , 5],
                         upper_quant = rw1_IV$summary.random$space.time[ , 6])
  
  rw2_I.df <- data.frame(x_axis = x_axis,
                         lower_quant = rw2_I$summary.random$space.time[ , 4], 
                         median = rw2_I$summary.random$space.time[ , 5],
                         upper_quant = rw2_I$summary.random$space.time[ , 6])
  
  rw2_II.df <- data.frame(x_axis = x_axis,
                          lower_quant = rw2_II$summary.random$space.time[ , 4], 
                          median = rw2_II$summary.random$space.time[ , 5],
                          upper_quant = rw2_II$summary.random$space.time[ , 6])
  
  rw2_III.df <- data.frame(x_axis = x_axis,
                           lower_quant = rw2_III$summary.random$space.time[ , 4], 
                           median = rw2_III$summary.random$space.time[ , 5],
                           upper_quant = rw2_III$summary.random$space.time[ , 6])
  
  rw2_IV.df <- data.frame(x_axis = x_axis,
                          lower_quant = rw2_IV$summary.random$space.time[ , 4], 
                          median = rw2_IV$summary.random$space.time[ , 5],
                          upper_quant = rw2_IV$summary.random$space.time[ , 6])
  
  
  
  
  library(latex2exp)
  rw1_I_plot <- ggplot(data = rw1_I.df, aes(x = x_axis)) + 
                  geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                              fill = "lightgrey", alpha = 0.6) +
                  geom_line(aes(y = median, col = "Median")) +
                  xlab("") + ylab("Type I Interaction") + ggtitle(TeX('$\\alpha_{t}$: RW1')) +
                  theme_bw() + 
                  labs(col = NULL) 
  
  rw2_I_plot <- ggplot(data = rw2_I.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("") + ylab("") + ggtitle(TeX('$\\alpha_{t}$: RW2')) +
    theme_bw() + 
    labs(col = NULL) 
  
  rw1_II_plot <- ggplot(data = rw1_II.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("") + ylab("Type II Interaction") + 
    theme_bw() + 
    labs(col = NULL) 
  
  rw2_II_plot <- ggplot(data = rw2_II.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("") + ylab("") + 
    theme_bw() + 
    labs(col = NULL) 
  
  rw1_III_plot <- ggplot(data = rw1_III.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("") + ylab("Type III Interaction") + 
    theme_bw() + 
    labs(col = NULL) 
  
  rw2_III_plot <- ggplot(data = rw2_III.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("") + ylab("") + 
    theme_bw() + 
    labs(col = NULL) 
  
  rw1_IV_plot <- ggplot(data = rw1_IV.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("space.time") + ylab("Type IV Interaction") + 
    theme_bw() + 
    labs(col = NULL) 
  
  rw2_IV_plot <- ggplot(data = rw2_IV.df, aes(x = x_axis)) + 
    geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                fill = "lightgrey", alpha = 0.6) +
    geom_line(aes(y = median, col = "Median")) +
    xlab("space.time") + ylab("") + 
    theme_bw() + 
    labs(col = NULL) 
  
  ggarrange(rw1_I_plot, rw2_I_plot,
            rw1_II_plot, rw2_II_plot,
            rw1_III_plot, rw2_III_plot,
            rw1_IV_plot, rw2_IV_plot,
            ncol = 2, nrow = 4,
            common.legend = TRUE, legend = "right")
  
  
  #temp <- 1:length(typeI$summary.random$space.time[ , 4])
  
  #typeI.df <- data.frame(x_axis = c(temp, temp, temp),
  #                       quantiles = c(typeI$summary.random$space.time[ , 4],
  #                                     typeI$summary.random$space.time[ , 5],
  #                                     typeI$summary.random$space.time[ , 6]),
  #                       type = c(rep("0.025 quantile", length(temp)),
  #                                rep("median", length(temp)),
  #                                rep("0.975 quantile", length(temp))))
  
  #typeI.df <- typeI.df %>% mutate(color = case_when(
  #  type == "median" ~ "black",
  #  TRUE ~ "grey"
  #))
  
  
  #ggplot(data = typeI.df) + 
  #  theme_bw() + 
  #  labs(col = NULL) +
  #  geom_line(aes(x = x_axis, y = quantiles, group = type, color = color)) +
  #  xlab("space.time") + ylab("Type I Interaction")
    
  
  
  #par(mfrow = c(2, 2))
  #matplot(typeI$summary.random$space.time[ , 4:6],
  #        lty=c(2,1,2), type="l", col=1,
  #        xlab = "space.time", ylab = "Type I Interaction")
  
  #matplot(typeII$summary.random$space.time[ , 4:6],
  #        lty=c(2,1,2), type="l", col=1,
  #        xlab = "space.time", ylab = "Type II Interaction effect")
  
  #matplot(typeIII$summary.random$space.time[ , 4:6],
  #        lty=c(2,1,2), type="l", col=1,
  #        xlab = "space.time", ylab = "Type III Interaction effect")
  
  #matplot(typeIV$summary.random$space.time[ , 4:6],
  #        lty=c(2,1,2), type="l", col=1,
  #        xlab = "space.time", ylab = "Type IV Interaction effect")
}




plot_proper_interaction <- function(proper_interaction, proper_full){
  #want to restructure to get year over county
  x_axis = 1:length(proper_interaction$summary.random$county[ , 5])
  
  #maybe sort the interactions
  #proper_interaction$summary.random$county = sort_proper_fitted(proper_interaction$summary.random$county, n, T)
  
  
  #Format for ggplot
  only_interaction.df <- data.frame(x_axis = x_axis,
                                 lower_quant = proper_interaction$summary.random$county[ , 4],
                                 median = proper_interaction$summary.random$county[ , 5],
                                 upper_quant = proper_interaction$summary.random$county[ , 6])
  
  
  interaction_and_more.df <- data.frame(x_axis = x_axis,
                                lower_quant = proper_full$summary.random$county.copy[ , 4],
                                median = proper_full$summary.random$county.copy[ , 5],
                                upper_quant = proper_full$summary.random$county.copy[ , 6])
  
  
  
  only_interaction_plot <- ggplot(data = only_interaction.df, aes(x = x_axis)) + 
                        geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                                    fill = "lightgrey", alpha = 0.6) +
                        geom_line(aes(y = median, col = "Median")) +
                        xlab("space.time") + ylab("Proper interaction") + 
                        theme_bw() + 
                        labs(col = NULL)
  
  interaction_and_more_plot <- ggplot(data = interaction_and_more.df, aes(x = x_axis)) + 
                        geom_ribbon(aes(ymin = lower_quant, ymax = upper_quant, col = "95% CI"),
                                    fill = "lightgrey", alpha = 0.6) +
                        geom_line(aes(y = median, col = "Median")) +
                        xlab("space.time") + ylab("Proper interaction") + 
                        theme_bw() + 
                        labs(col = NULL)
  
  ggarrange(only_interaction_plot, interaction_and_more_plot,
            ncol = 1, nrow = 2,
            common.legend = TRUE, legend = "right")
  
  #par(mfrow = c(2, 1))
  #matplot(proper_interaction$summary.random$county[ , 5],
  #                lty=c(1), type="l", col=1,
  #                xlab = "", ylab = "proper interaction")
  
  #matplot(proper_full$summary.random$county.copy[ , 5],
  #        lty=c(1), type="l", col=1,
  #        xlab = "", ylab = "proper interaction")
  
}




#Plot precision of interaction RW1
plot_std_interactions_RW1 <- function(rw1_typeI, rw1_typeII, 
                                  rw1_typeIII, rw1_typeIV){
  
  #Transform from precision to standard deviation
  rw1_std_I <- inla.tmarginal(function(x) sqrt(1/x),
                              rw1_typeI$marginals.hyperpar$`Precision for space.time`)
  
  rw1_std_II <- inla.tmarginal(function(x) sqrt(1/x),
                               rw1_typeII$marginals.hyperpar$`Precision for space.time`)
  
  rw1_std_III <- inla.tmarginal(function(x) sqrt(1/x),
                                rw1_typeIII$marginals.hyperpar$`Precision for space.time`)
  
  rw1_std_IV <- inla.tmarginal(function(x) sqrt(1/x),
                               rw1_typeIV$marginals.hyperpar$`Precision for space.time`)
  
  
  #Format for ggplot
  rw1_std_df <- data.frame(x_axis = c(rw1_std_I[, 1], rw1_std_II[, 1],
                                      rw1_std_III[, 1], rw1_std_IV[, 1]),
                           y_axis = c(rw1_std_I[, 2], rw1_std_II[, 2],
                                      rw1_std_III[, 2], rw1_std_IV[, 2]),
                           type = c(rep("Type I", length(rw1_std_I[, 1])),
                                    rep("Type II", length(rw1_std_II[, 1])),
                                    rep("Type III", length(rw1_std_III[, 1])),
                                    rep("Type IV", length(rw1_std_IV[, 1]))))
  
  #Plot
  gplt <- ggplot(data=rw1_std_df,
                aes(x=x_axis, group=type, fill=type)) +
                geom_density(adjust=1.5, alpha=.8) +
                theme_bw() + 
                theme(axis.title=element_text(size=14)) +
                labs(fill = NULL) +
                xlab(expression(sigma)) + ylab(expression(f(sigma)))
  plt <- ggarrange(gplt, common.legend = TRUE, legend = "top")
  annotate_figure(plt, top = text_grob("Interaction hyperparameters for Improper_1 models", 
                                         color = "black", size = 14))
  
}

#Plot precision of interaction RW2
plot_std_interactions_RW2 <- function(typeI,
                                      typeII, 
                                      typeIII, 
                                      typeIV,
                                      title){
  
  #Transform from precision to standard deviation
  std_I <- inla.tmarginal(function(x) sqrt(1/x),
                          typeI$marginals.hyperpar$`Precision for space.time`)
  
  std_II <- inla.tmarginal(function(x) sqrt(1/x),
                           typeII$marginals.hyperpar$`Precision for space.time`)
  
  std_III <- inla.tmarginal(function(x) sqrt(1/x),
                            typeIII$marginals.hyperpar$`Precision for space.time`)
  
  std_IV <- inla.tmarginal(function(x) sqrt(1/x),
                           typeIV$marginals.hyperpar$`Precision for space.time`)
  
  
  #Format for ggplot
  std_df_I_III <- data.frame(x_axis = c(std_I[, 1], std_III[, 1]),
                             y_axis = c(std_I[, 2], std_III[, 2]),
                             type = c(rep("Type I", length(std_I[, 1])),
                                      rep("Type III", length(std_III[, 1]))))
  
  std_df_II_IV <- data.frame(x_axis = c(std_II[, 1], std_IV[, 1]),
                             y_axis = c(std_II[, 2], std_IV[, 2]),
                             type = c(rep("Type II", length(std_II[, 1])),
                                      rep("Type IV", length(std_IV[, 1]))))
  
  #Plot
  I_III <- ggplot(data=std_df_I_III,
                  aes(x=x_axis, group=type, fill=type)) +
                  geom_density(adjust=1.5, alpha=.8) +
                  theme_bw() +
                  labs(fill = NULL) +
                  theme(axis.title=element_text(size=14)) +
                  xlab(expression(sigma)) + ylab(expression(f(sigma)))
  
  II_IV <- ggplot(data=std_df_II_IV,
                  aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.8) +
    theme_bw() +
    labs(fill = NULL) +
    theme(axis.title=element_text(size=14)) +
    xlab(expression(sigma)) + ylab(expression(f(sigma)))
  
  plt <- ggarrange(I_III, II_IV, ncol = 2, nrow = 1, 
            common.legend = FALSE, legend = "top")
  
  annotate_figure(plt, top = text_grob("Interaction hyperparameters for Improper_2 models", 
                                       color = "black", size = 12))
  
}




#Plot precision of proper interaction
plot_proper_hyperparameters <- function(interaction_only, full){
  #Transform from precision to standard deviation
  std_interaction_only <- inla.tmarginal(function(x) sqrt(1/x),
                                 interaction_only$marginals.hyperpar$`Precision for county`)
  
  std_interaction_full <- inla.tmarginal(function(x) sqrt(1/x),
                                         full$marginals.hyperpar$`Precision for county.copy`)
  
  std.df <- data.frame(x_axis = c(std_interaction_only[, 1], std_interaction_full[, 1]),
                       y_axis = c(std_interaction_only[, 2], std_interaction_full[, 2]),
                       type = c(rep("Proper_onlyInt", length(std_interaction_only[, 1])),
                                rep("Proper_full", length(std_interaction_full[, 1]))))
  
  std_plot <- ggplot(data=std.df,
                     aes(x=x_axis, group=type, fill=type)) +
                     geom_density(adjust=1.5, alpha=.8) +
                     theme_bw() +
                     labs(fill = NULL) +
                     theme(axis.title=element_text(size=14)) +
                     xlab(expression(sigma)) + ylab(TeX("$f(*)$"))
  
  
  interaction_lambda <- interaction_only$marginals.hyperpar$`Lambda for county`
  full_lambda <- full$marginals.hyperpar$`Lambda for county.copy`
  
  lambda.df <- data.frame(x_axis = c(interaction_lambda[, 1], full_lambda[, 1]),
                          y_axis = c(interaction_lambda[, 2], full_lambda[, 2]),
                          type = c(rep("Proper_onlyInt", length(interaction_lambda[, 1])),
                                   rep("Proper_full", length(full_lambda[, 1]))))
  
  lambda_plot <- ggplot(data=lambda.df,
                        aes(x=x_axis, group=type, fill=type)) +
                        geom_density(adjust=3, alpha=.8) +
                        theme_bw() +
                        labs(fill = NULL) +
                        theme(axis.title=element_text(size=14)) +
                        xlab(expression(lambda)) + ylab(NULL)
  
  
  interaction_rho <- interaction_only$marginals.hyperpar$`GroupRho for county`
  full_rho <- full$marginals.hyperpar$`GroupRho for county.copy`
  
  rho.df <- data.frame(x_axis = c(interaction_rho[, 1], full_rho[, 1]),
                       y_axis = c(interaction_rho[, 2], full_rho[, 2]),
                       type = c(rep("Proper_onlyInt", length(interaction_rho[, 1])),
                                rep("Proper_full", length(full_rho[, 1]))))
  
  
  rho_plot <- ggplot(data=rho.df,
                        aes(x=x_axis, group=type, fill=type)) +
                        geom_density(adjust=5, alpha=.8) +
                        theme_bw() +
                        labs(fill = NULL) +
                        theme(axis.title=element_text(size=14)) +
                        xlab(expression(rho)) + ylab(NULL)
  
  
  plt <- ggarrange(std_plot, lambda_plot, rho_plot,
                    ncol = 3, nrow = 1,
                    common.legend = TRUE, legend = "top")
  
  annotate_figure(plt, top = text_grob("Interaction hyperparameters for Proper models", 
                                       color = "black", size = 15))
}







#Plot fitted values against actual values (all together)
plot_fitted_vs_actual_together <- function(actual_data,
                                           base_model,
                                           typeII_model,
                                           n, T){
  
  #Format for ggplot
  fitted.df <- data.frame(lower_quant = typeII_model$summary.fitted.values$'0.025quant',
                          upper_quant = typeII_model$summary.fitted.values$'0.975quant',
                          x_axis = 1:(n*T),
                          actual = actual_data$rate)
  
  base.df <- data.frame(lower_quant = base_model$summary.fitted.values$'0.025quant',
                          upper_quant = base_model$summary.fitted.values$'0.975quant',
                          x_axis = 1:(n*T),
                          actual = actual_data$rate)
  
  
  II <- ggplot(data = fitted.df) + 
    geom_ribbon(data = fitted.df, 
                aes(x = x_axis,
                    ymin = lower_quant,
                    ymax = upper_quant, col = "95% CI"),
                fill = "black",
                alpha = 0.2) +
    geom_point(aes(x = x_axis, y = actual, col = "True values")) + 
    xlab("Datum ID") +ylab("Relative Risk") + 
    ggtitle("Model with type II interaction\n Fitted Relative Risk against True Rate") +
    theme_bw()
  
  base <- ggplot(data = base.df) + 
    geom_ribbon(aes(x = x_axis, 
                    ymin = lower_quant,
                    ymax = upper_quant, 
                    col = "95% CI"),
                fill = "black") +
    geom_point(aes(x = x_axis, y = actual, col = "True values")) + 
    xlab("Datum ID") +ylab("Relative Risk") + 
    ggtitle("Model without interactions\n Fitted Relative Risk against True Rate") +
    theme_bw()
  
  ggarrange(base, II, ncol = 2, common.legend = TRUE, legend = "right")
  
}



#Plot a time series for a county showing fitted values vs actual values
county_time_series <- function(actual,
                               fitted_model,
                               county,
                               n, T,
                               title = TRUE,
                               xlab,
                               ylab){
  years <- 1968:1988
  values.df <- data.frame(years = years, 
                          true_rate = actual$rate[seq(county, n*T, by = n)] * 1E5,
                          fitted_rate = fitted_model$summary.fitted.values[seq(county,n*T,by=n), 4] * 1E5,
                          lower_quant = fitted_model$summary.fitted.values[seq(county,n*T,by=n), 3] * 1E5,
                          upper_quant = fitted_model$summary.fitted.values[seq(county,n*T,by=n), 5] * 1E5)
  if(title){
  plt <- ggplot(data = values.df, aes(x = years)) + 
          geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                      fill = "#F8766D", alpha = 0.6) +
          geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
          geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
          xlab(xlab) + ylab(ylab) + ggtitle(actual[county, ]$name) +
          labs(col = NULL) +
          theme_bw() + 
          theme(axis.title=element_text(size=14))
  } else {
    plt <- ggplot(data = values.df, aes(x = years)) + 
      geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                  fill = "#F8766D", alpha = 0.6) +
      geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
      geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
      xlab(xlab) + ylab(ylab) +
      labs(col = NULL) +
      theme_bw() +
      theme(axis.title=element_text(size=14))
  }
  plt <- plt + scale_color_manual(values=c("#F8766D", "black", "#00BFC4"))
  return(plt)
}

select_county_timeseries <- function(actual_data,
                                     base_model, 
                                     RW1_II,
                                     proper_interaction,
                                     counties,
                                     n, T){
  base_plt1 <- county_time_series(actual_data, base_model, counties[1], n, T, title = T, xlab = "", ylab = "Rate pr. 100000")
  base_plt2 <- county_time_series(actual_data, base_model, counties[2], n, T, title = T, xlab = "", ylab = "")
  base_plt3 <- county_time_series(actual_data, base_model, counties[3], n, T, title = T, xlab = "", ylab = "")
  base_plt4 <- county_time_series(actual_data, base_model, counties[4], n, T, title = T, xlab = "", ylab = "")
  
  II_plt1 <- county_time_series(actual_data, RW1_II, counties[1], n, T, title = F, xlab = "", ylab = "Rate pr. 100000")
  II_plt2 <- county_time_series(actual_data, RW1_II, counties[2], n, T, title = F, xlab = "", ylab = "")
  II_plt3 <- county_time_series(actual_data, RW1_II, counties[3], n, T, title = F, xlab = "", ylab = "")
  II_plt4 <- county_time_series(actual_data, RW1_II, counties[4], n, T, title = F, xlab = "", ylab = "")
  
  #Proper models must be sorted first
  proper_interaction$summary.fitted.values <- sort_proper_fitted(proper_interaction$summary.fitted.values,
                                                                 n, T)
  
  proper_plt1 <- county_time_series(actual_data, proper_interaction, counties[1], n, T, title = F, xlab = "year", ylab = "Rate pr. 100000")
  proper_plt2 <- county_time_series(actual_data, proper_interaction, counties[2], n, T, title = F, xlab = "year", ylab = "")
  proper_plt3 <- county_time_series(actual_data, proper_interaction, counties[3], n, T, title = F, xlab = "year", ylab = "")
  proper_plt4 <- county_time_series(actual_data, proper_interaction, counties[4], n, T, title = F, xlab = "year", ylab = "")
  
  ggarrange(base_plt1, base_plt2, base_plt3, base_plt4,
            II_plt1, II_plt2, II_plt3, II_plt4,
            proper_plt1, proper_plt2, proper_plt3, proper_plt4,
            ncol = 4, nrow = 3, 
            common.legend = TRUE, legend = "top")
  
}








#Plot heatmap for 1 year
case_count_plot_1_year <- function(sf_data, 
                                   year,
                                   hardcoded_bins,
                                   title){
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
  #scale = heat.colors(8, rev= TRUE)
  
  p <- ggplot(data = sf_data[sf_data$year == year, ]) + 
      geom_sf(aes(fill = rate), #Plots death rate
              alpha = 1,
              color="black") + 
      ggtitle(title) +
    theme(plot.title = element_text(size = 14),
            axis.title.x = element_blank(), #Remove axis and background grid
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.background = element_blank(),
            panel.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
            plot.margin =  unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
            legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "cm"),
            panel.spacing = unit(1, 'lines')) +
      guides(fill=guide_legend(title=NULL, reverse = TRUE, label.position = "right")) + #Remove colorbar title
      binned_scale( #Scaling the color
        aesthetics = "fill",
        scale_name = "gradientn",
        palette = function(x) c(scale),
        labels = function(x){x},
        breaks = hardcoded_bins,
        limits = c(0, 100),
        guide = "colorscale")
  
  
  return(p)
}







violin_plot_rate <- function(rates.df){
  ggplot(data = rates.df, aes(x=type, y=fitted_rates, fill=type)) + 
  geom_violin() +
  xlab("Model") +
  theme(legend.position="none") +
    theme_bw() + 
  ylab("Relative Risk")
}


####
extract_res_inla = function(result, #Fitted model, i suppose?
                            n_samples = 10000,
                            return_hist_data = T,
                            check_samples = 100){
  
  #n_fixed = ifelse(is.na(nrow(result$summary.fixed)), 0, nrow(result$summary.fixed))
  #stopifnot(n_fixed > 0)
  
  posterior_samples = inla.hyperpar.sample(n_samples, result) #Sample the posterior
  
  
  #effect_names = word(colnames(posterior_samples), start = -1) # Extract names
  #stopifnot(length(effect_names) > 0) # Need random effect to work
  
  
  rho_samples = data.frame(matrix(ncol = length(effect_names), nrow = n_samples))
  
  
  colnames(rho_samples) = effect_names
  
  #Transform precision to variance
  for (i in 1:n_samples) {
    rho_samples[i,] = as.numeric(1/posterior_samples[i,])
  }
  sigma = sqrt(rowSums(rho_samples)) # Square root total marginal variance
  rho_samples = rho_samples/rowSums(rho_samples) #Calculate proportion of variance in each effect
  
  #Create dataframe to store the results
  inference_results = data.frame(matrix(ncol = 3, nrow = 1 + n_fixed + length(effect_names)))
  colnames(inference_results) = c("median", "0.025quant", "0.975quant")
  
  #These first rows are the fixed effects
  for (i in 1:n_fixed) {
    inference_results[i, ] = result$summary.fixed[c("0.5quant", "0.025quant", "0.975quant")][i,]
  }
  
  rownames(inference_results)[1:n_fixed] = rownames(result$summary.fixed)
  
  for (i in 1:length(effect_names)) {
    effect = effect_names[i]
    cred_int = ci(rho_samples[effect], method = "ETI")
    inference_results[n_fixed + i,] = c(median(rho_samples[effect][,1]), cred_int$CI_low, cred_int$CI_high)
  }
  
  rownames(inference_results)[(n_fixed + 1):(n_fixed + length(effect_names))] = paste0("rho_", effect_names)
  
  
  #Create equi-tailed CI
  cred_int_sigma = ci(sigma, method ="ETI")
  inference_results[n_fixed + length(effect_names) + 1,] = c(median(sigma), cred_int_sigma$CI_low, cred_int_sigma$CI_high)
  rownames(inference_results)[n_fixed + length(effect_names) + 1] = "sigma"
  
  # Include the samples in the returned object. For making posterior densities
  if (return_hist_data){
    hist_data = data.frame(matrix(ncol = length(effect_names) + 1, nrow = n_samples))
    colnames(hist_data)[1:length(effect_names)] = effect_names
    colnames(hist_data)[length(effect_names)+1] = "sigma"
    for (i in 1:length(effect_names)) {
      hist_data[effect_names[i]][,1] = rho_samples[,i]
    }
    colnames(hist_data)[1:length(effect_names)] = paste0("rho_", effect_names)
    hist_data$sigma = as.numeric(sigma)
    inference_results = list(inference_results = inference_results, hist_data = hist_data)
  }
  
  # Calculate posterior probabilities 
  inference_results$posterior_prob = c(paste("Posterior probability of", effect_names[1], ">",effect_names[2], sep =" "), mean(ifelse(rho_samples[,1] > rho_samples[,2], 1, 0)))
  
  return(inference_results)
}



