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

abs_error <- function(x, mu = 0){
  return(abs(x - mu))
}

find_CRPS_ae_one_year <- function(marginals, pop, true_values, year){
  crps = rep(0, n); abs_error = rep(0, n)
  mean_marg = rep(0, n); mean_predicted = rep(0, n)
  var_marg = rep(0, n); var_predicted = rep(0, n)
  
  for(i in 1:n){
    mean_marg[i] = find_mean_marginal(marginals[year, i])
    var_marg[i] = find_var_marginal(marginals[year, i])
    
    mean_predicted[i] = pop[((year - 1) * n + i)] * mean_marg[i]
    var_predicted[i] = pop[((year - 1) * n + i)] * mean_marg[i] + 
      (pop[((year - 1) * n + i)]**2) * var_marg[i] 
    
    crps[i] = crpsNormal(true_values[((year - 1) * n + i)], 
                         mu = mean_predicted[i],
                         sig = sqrt(var_predicted[i]))
    
    abs_error[i] = abs_error(true_values[((year - 1) * n + i)],
                             mu = mean_predicted[i])
  }
  return(list(crps = crps, 
              ae = abs_error,
              mean_marg = mean_marg,
              var_marg = var_marg,
              mean_pred = mean_predicted,
              var_pred = var_predicted))
}


find_CRPS_ae_all_years <- function(marginals, pop, true){
  avg_crps_each_year = rep(0, 10)
  avg_ae_each_year = rep(0, 10)
  for(t in 1:10){
    temp                  = find_CRPS_ae_one_year(marginals, pop, true, t)
    avg_crps_each_year[t] = mean(temp$crps)
    avg_ae_each_year[t] = mean(temp$ae)
  }
  return(list(yearly_avg_crps = avg_crps_each_year,
              yearly_avg_ae = avg_ae_each_year))
}



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
              ggtitle("Random effect")
  
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
                    ggtitle("Random and fixed effect")
    
    
  ggarrange(fixed_plot, ar1_plot, ar1_fixed_plot, 
            ncol = 3, nrow = 1)
}



#Plot posterior intercept, temporal effects and spatial effects
plot_spatial_effects <- function(improper,
                                 proper,
                                 map,
                                 n){
  #Function that produces four plots: The posterior intercept, 
  #posterior structured temporal effect along with 2.5% and 97.5% quantiles
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale_1 = scale_col[c(3,10,13,17,21,24,27,30)] #Select color scale to be more red
  
  improper_mean <- improper$summary.random$county$mean[(n+1):(2*n)]
  proper_mean <- proper$summary.random$county$mean
  
  temp_ohio_map <- map[ ,c("geometry", "NAME")]
  temp_ohio_map$improper_mean <- improper_mean
  temp_ohio_map$proper_mean <- proper_mean
  
  #Set the theme to minimal
  theme_set(theme(panel.background = element_blank()))
  improper_mean_plot <- ggplot(data = temp_ohio_map) + 
                        geom_sf(aes(fill = improper_mean), 
                                alpha = 1,
                                color="black") + ggtitle("Mean Besag") +
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
  
  proper_mean_plot <- ggplot(data = temp_ohio_map) + 
                      geom_sf(aes(fill = improper_mean), 
                              alpha = 1,
                              color="black") + ggtitle("Mean proper Besag") +
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
  scale_1 = scale_col[c(3,10,13,17,21,24,27,30)] #Select color scale to be more red
  
  improper_sd <- improper$summary.random$county$sd[(n+1):(2*n)]
  proper_sd <- proper$summary.random$county$sd
  
  temp_ohio_map <- map[ ,c("geometry", "NAME")]
  temp_ohio_map$improper_sd <- improper_sd
  temp_ohio_map$proper_sd <- proper_sd
  
  #Set the theme to minimal
  theme_set(theme(panel.background = element_blank()))
  improper_sd_plot <- ggplot(data = temp_ohio_map) + 
    geom_sf(aes(fill = improper_sd), 
            alpha = 1,
            color="black") + ggtitle("Standard deviation Besag") +
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
  
  proper_sd_plot <- ggplot(data = temp_ohio_map) + 
    geom_sf(aes(fill = proper_sd), 
            alpha = 1,
            color="black") + ggtitle("Standard deviation proper Besag") +
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

  
  ggarrange(improper_sd_plot, proper_sd_plot,
            ncol = 2, nrow = 1,
            common.legend = FALSE)
  
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
    xlab(expression(sigma)) + ylab(expression(f(sigma))) 
  
  phi_temporal_plot <- ggplot(data=phi_temporal_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=6.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(lambda)) + ylab(expression(f(lambda)))
  
  ggarrange(std_temporal_plot, phi_temporal_plot,
            common.legend = T, legend = "right")
  
  
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
    xlab(expression(sigma)) + ylab(expression(f(sigma)))
  
  phi_spatial_plot <- ggplot(data=phi_spatial_df,
                             aes(x=x_axis)) +
    stat_density(fill = "#00BFC4", adjust=3, alpha=.8,
                 position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(lambda)) + ylab(expression(f(lambda)))
  
  ggarrange(std_spatial_plot, phi_spatial_plot)
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
                                type = c(rep("Base", length(std_ar1_base[, 1])),
                                         rep("Full", length(std_ar1_full[, 1]))))
  
  
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
                 adjust=1.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) 
  
  rho_temporal_plot <- ggplot(data=rho_temporal_df) + 
    stat_density(aes(x=x_axis, group=type, fill=type),
                 adjust=6.5, alpha=.8, position = "identity") +
    theme_bw() +
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(rho)) + ylab(expression(f(rho)))
  
  ggarrange(std_temporal_plot, rho_temporal_plot,
            common.legend = T, legend = "right")
  
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
                               type = c(rep("Base", length(std_CAR_base[, 1])),
                                        rep("Full", length(std_CAR_full[, 1]))))
  
  
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
  
  ggarrange(std_spatial_plot, lambda_spatial_plot,
            common.legend = T, legend = "right")
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
  ggplot(data=rw1_std_df,
         aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.8) +
    theme_bw() + 
    theme(axis.title=element_text(size=14)) +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(expression(f(sigma)))
  
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
  
  ggarrange(I_III, II_IV, ncol = 2, nrow = 1, 
            common.legend = FALSE)
  
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
                       type = c(rep("Interaction only", length(std_interaction_only[, 1])),
                                rep("Full model", length(std_interaction_full[, 1]))))
  
  std_plot <- ggplot(data=std.df,
                     aes(x=x_axis, group=type, fill=type)) +
                     geom_density(adjust=1.5, alpha=.8) +
                     theme_bw() +
                     labs(fill = NULL) +
                     theme(axis.title=element_text(size=14)) +
                     xlab(expression(sigma)) + ylab(expression(f(sigma)))
  
  
  interaction_lambda <- interaction_only$marginals.hyperpar$`Lambda for county`
  full_lambda <- full$marginals.hyperpar$`Lambda for county.copy`
  
  lambda.df <- data.frame(x_axis = c(interaction_lambda[, 1], full_lambda[, 1]),
                          y_axis = c(interaction_lambda[, 2], full_lambda[, 2]),
                          type = c(rep("Interaction only", length(interaction_lambda[, 1])),
                                   rep("Full model", length(full_lambda[, 1]))))
  
  lambda_plot <- ggplot(data=lambda.df,
                        aes(x=x_axis, group=type, fill=type)) +
                        geom_density(adjust=3, alpha=.8) +
                        theme_bw() +
                        labs(fill = NULL) +
                        theme(axis.title=element_text(size=14)) +
                        xlab(expression(lambda)) + ylab(expression(f(lambda)))
  
  
  interaction_rho <- interaction_only$marginals.hyperpar$`GroupRho for county`
  full_rho <- full$marginals.hyperpar$`GroupRho for county.copy`
  
  rho.df <- data.frame(x_axis = c(interaction_rho[, 1], full_rho[, 1]),
                       y_axis = c(interaction_rho[, 2], full_rho[, 2]),
                       type = c(rep("Interaction only", length(interaction_rho[, 1])),
                                rep("Full model", length(full_rho[, 1]))))
  
  
  rho_plot <- ggplot(data=rho.df,
                        aes(x=x_axis, group=type, fill=type)) +
                        geom_density(adjust=5, alpha=.8) +
                        theme_bw() +
                        labs(fill = NULL) +
                        theme(axis.title=element_text(size=14)) +
                        xlab(expression(rho)) + ylab(expression(f(rho)))
  
  
  ggarrange(std_plot, lambda_plot, rho_plot,
            ncol = 3, nrow = 1,
            common.legend = TRUE, legend = "right")
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
          theme_bw()
  } else {
    plt <- ggplot(data = values.df, aes(x = years)) + 
      geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                  fill = "#F8766D", alpha = 0.6) +
      geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
      geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
      xlab(xlab) + ylab(ylab) +
      labs(col = NULL) +
      theme_bw()
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
            common.legend = TRUE, legend = "right")
  
}








#Plot heatmap for 1 year
case_count_plot_1_year <- function(sf_data, 
                                   year,
                                   hardcoded_bins){
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
  #scale = heat.colors(8, rev= TRUE)
  
  p <- ggplot(data = sf_data[sf_data$year == year, ]) + 
      geom_sf(aes(fill = rate), #Plots death rate
              alpha = 1,
              color="black") + 
      ggtitle(year) +
      theme(plot.title = element_text(size = 10),
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






