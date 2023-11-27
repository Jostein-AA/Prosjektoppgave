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
      A <- t(eigens$vectors[ ,(nrow(eigens$vectors) - 2 * n - T + 3):nrow(eigens$vectors)])
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


#Plot the intercept
plot_intercept <- function(fitted_model){
  #Format intercept for ggplot
  intercept.df <- data.frame(x = fitted_model$marginals.fixed$`(Intercept)`[, 1])
  intercept.df$y <- fitted_model$marginals.fixed$`(Intercept)`[, 2]
  
  ggplot(data = intercept.df, aes(x = x, y = y)) + 
    theme_bw() +
    geom_line() + 
    xlab(expression(mu)) + ylab(expression(f(mu))) + 
    ggtitle("Posterior Density of Intercept")
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
    geom_line() + 
    geom_line(data = temporal_RW1.df, aes(years, lower_quant), linetype = "dashed") + 
    geom_line(data = temporal_RW1.df, aes(years, upper_quant), linetype = "dashed") + 
    xlab("year: t") + ylab(expression(alpha[t])) + 
    ggtitle("Structured Temporal Effect RW1")
  
  temporal_RW2 <- ggplot(data = temporal_RW2.df, aes(years, median)) + 
    theme_bw() +
    geom_line() + 
    geom_line(data = temporal_RW2.df, aes(years, lower_quant), linetype = "dashed") + 
    geom_line(data = temporal_RW2.df, aes(years, upper_quant), linetype = "dashed") + 
    xlab("year: t") + ylab(expression(alpha[t])) + 
    ggtitle("Structured Temporal Effect RW2")
  
  ggarrange(temporal_RW1, temporal_RW2)
}


#Plot posterior intercept, temporal effects and spatial effects
plot_spatial_effects <- function(fitted_model, map, n){
  #Function that produces four plots: The posterior intercept, 
  #posterior structured temporal effect along with 2.5% and 97.5% quantiles
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale_1 = scale_col[c(3,10,13,17,21,24,27,30)] #Select color scale to be more red
  
  spatial_structured_effect_mean <- fitted_model$summary.random$county$mean[(n+1):(2*n)]
  spatial_structured_effect_sd <- fitted_model$summary.random$county$sd[(n+1):(2*n)]
  temp_ohio_map <- map[ ,c("geometry", "NAME")]
  temp_ohio_map$mean <- spatial_structured_effect_mean
  temp_ohio_map$sd <- spatial_structured_effect_sd
  
  #Set the theme to minimal
  theme_set(theme(panel.background = element_blank()))
  p_1 <- ggplot(data = temp_ohio_map) + 
    geom_sf(aes(fill = mean), 
            alpha = 1,
            color="black") + ggtitle("Mean Spatial Structured Effect each County") +
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
  
  p_2 <- ggplot(data = temp_ohio_map) + 
    geom_sf(aes(fill = sd), 
            alpha = 1,
            color="black") + ggtitle("Standard deviation of Structured Spatial Effect\n for each county") +
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
  
  ggarrange(p_1, p_2,
            ncol = 2, nrow = 1,
            common.legend = FALSE)
}


#Plot posterior hyperparameters of temporal and random effects
plot_temporal_spatial_hyperparameters <- function(fitted_RW1, fitted_RW2){
  #inla.tmarginal to transform from precision to standard deviation
  std_year_RW1 <- inla.tmarginal(function(x) sqrt(1/x),
                                 fitted_RW1$marginals.hyperpar$`Precision for year`)
  
  std_year_RW2 <- inla.tmarginal(function(x) sqrt(1/x),
                                 fitted_RW2$marginals.hyperpar$`Precision for year`)
  
  std_county_RW1 <- inla.tmarginal(function(x) sqrt(1/x),
                                   fitted_RW1$marginals.hyperpar$`Precision for county`)
  
  std_county_RW2 <- inla.tmarginal(function(x) sqrt(1/x),
                                   fitted_RW2$marginals.hyperpar$`Precision for county`)
  
  
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
                                             length(fitted_RW2$marginals.hyperpar$`Phi for year`[, 1])))
  )
  
  
  std_spatial_df <- data.frame(x_axis = c(std_county_RW1[, 1], std_county_RW2[, 1]),
                                y_axis = c(std_county_RW1[, 2], std_county_RW2[, 2]),
                                type = c(rep("RW1", length(std_county_RW1[, 1])),
                                         rep("RW2", length(std_county_RW2[, 1]))))
  
  phi_spatial_df <- data.frame(x_axis = c(fitted_RW1$marginals.hyperpar$`Phi for county`[, 1],
                                           fitted_RW2$marginals.hyperpar$`Phi for county`[, 1]),
                                y_axis = c(fitted_RW1$marginals.hyperpar$`Phi for county`[, 2],
                                           fitted_RW2$marginals.hyperpar$`Phi for county`[, 2]),
                                type = c(rep("RW1", 
                                             length(fitted_RW1$marginals.hyperpar$`Phi for county`[, 1])),
                                         rep("RW2", 
                                             length(fitted_RW2$marginals.hyperpar$`Phi for county`[, 1])))
  )

  
  # With transparency (right)
  std_temporal_plot <- ggplot(data=std_temporal_df,
                              aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) + 
    ggtitle("Posterior Density of Standard deviation\n of Temporal Random Effects") 
  
  phi_temporal_plot <- ggplot(data=phi_temporal_df,
                              aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    labs(fill = NULL) +
    xlab(expression(phi)) + ylab(expression(f(phi))) + 
    ggtitle("Posterior Density of Temporal Mixing Parameter") 
  
  std_spatial_plot <- ggplot(data=std_spatial_df,
                              aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    labs(fill = NULL) +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) + 
    ggtitle("Posterior Density of Standard deviation\n of Spatial Random Effects")
  
  phi_spatial_plot <- ggplot(data=phi_spatial_df,
                              aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    labs(fill = NULL) +
    xlab(expression(phi)) + ylab(expression(f(phi))) + 
    ggtitle("Posterior Density of Spatial Mixing Parameter") 
  
  
  ggarrange(std_temporal_plot, phi_temporal_plot,
            std_spatial_plot, phi_spatial_plot,
            common.legend = T, legend = "right")
  
  
}



#Plot the interactions
plot_interaction <- function(typeI, typeII, typeIII, typeIV){
  par(mfrow = c(2, 2))
  matplot(typeI$summary.random$space.time[ , 4:6],
          lty=c(2,1,2), type="l", col=1,
          xlab = "Datum ID", ylab = "Type I Interaction effect")
  
  matplot(typeII$summary.random$space.time[ , 4:6],
          lty=c(2,1,2), type="l", col=1,
          xlab = "Datum ID", ylab = "Type II Interaction effect")
  
  matplot(typeIII$summary.random$space.time[ , 4:6],
          lty=c(2,1,2), type="l", col=1,
          xlab = "Datum ID", ylab = "Type III Interaction effect")
  
  matplot(typeIV$summary.random$space.time[ , 4:6],
          lty=c(2,1,2), type="l", col=1,
          xlab = "Datum ID", ylab = "Type IV Interaction effect")
}



#Plot precision of interaction RW1
plot_std_interactions_RW1 <- function(typeI,
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
  std_df <- data.frame(x_axis = c(std_I[, 1], std_II[, 1], std_III[, 1], std_IV[, 1]),
                       y_axis = c(std_I[, 2], std_II[, 2], std_III[, 2], std_IV[, 2]),
                       type = c(rep("Type I", length(std_I[, 1])),
                                rep("Type II", length(std_II[, 1])),
                                rep("Type III", length(std_III[, 1])),
                                rep("Type IV", length(std_IV[, 1]))))
  
  #Plot
  ggplot(data=std_df,
         aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) + 
    ggtitle("Posterior Density of Standard deviation\n of Interaction")
}

#Plot precision of interaction
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
                  geom_density(adjust=1.5, alpha=.4) +
                  theme_bw() +
                  xlab(expression(sigma)) + ylab(expression(f(sigma))) + 
                ggtitle("Posterior Density of Standard deviation\n of Interaction")
  
  II_IV <- ggplot(data=std_df_II_IV,
                  aes(x=x_axis, group=type, fill=type)) +
    geom_density(adjust=1.5, alpha=.4) +
    theme_bw() +
    xlab(expression(sigma)) + ylab(expression(f(sigma))) + 
    ggtitle("Posterior Density of Standard deviation\n of Interaction")
  
  ggarrange(I_III, II_IV, ncol = 2, nrow = 1, common.legend = FALSE)
  
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



#Plot a time series for every county showing fitted values vs actual values
every_county_time_series <- function(fitted_model,
                                     actual,
                                     counties,
                                     n, T){
  years <- 1968:1988
  #Format for plotting
  
  c = counties[1]
  county1 <- data.frame(years = years, 
                        true_rate = actual$rate[seq(c, n*T, by = n)],
                        fitted_rate = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 4],
                        lower_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 3],
                        upper_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 5])
  
  c = counties[2]
  county2 <- data.frame(years = years, 
                        true_rate = actual$rate[seq(c, n*T, by = n)],
                        fitted_rate = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 4],
                        lower_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 3],
                        upper_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 5])
  
  c = counties[3]
  county3 <- data.frame(years = years, 
                        true_rate = actual$rate[seq(c, n*T, by = n)],
                        fitted_rate = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 4],
                        lower_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 3],
                        upper_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 5])
  
  c = counties[4]
  county4 <- data.frame(years = years, 
                        true_rate = actual$rate[seq(c, n*T, by = n)],
                        fitted_rate = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 4],
                        lower_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 3],
                        upper_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 5])
  
  c = counties[5]
  county5 <- data.frame(years = years, 
                        true_rate = actual$rate[seq(c, n*T, by = n)],
                        fitted_rate = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 4],
                        lower_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 3],
                        upper_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 5])
  
  c = counties[6]
  county6 <- data.frame(years = years, 
                        true_rate = actual$rate[seq(c, n*T, by = n)],
                        fitted_rate = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 4],
                        lower_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 3],
                        upper_quant = fitted_model$summary.fitted.values[seq(c,n*T,by=n), 5])
  
 
  
  plt1 <- ggplot(data = county1, aes(x = years)) + 
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "pink", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
    geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
    xlab("year") + ylab("Relative Risk") + ggtitle(actual[counties[1], ]$name) +
    labs(col = NULL) +
    theme_bw()
  
  plt2 <- ggplot(data = county2, aes(x = years)) + 
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "pink", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
    geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
    xlab("year") + ylab("Relative Risk") + ggtitle(actual[counties[2], ]$name) +
    labs(col = NULL) +
    theme_bw()
  
  plt3 <- ggplot(data = county3, aes(x = years)) + 
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "pink", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
    geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
    xlab("year") + ylab("Relative Risk") + ggtitle(actual[counties[3], ]$name) +
    labs(col = NULL) +
    theme_bw()
  
  plt4 <- ggplot(data = county4, aes(x = years)) + 
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "pink", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
    geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
    xlab("year") + ylab("Relative Risk") + ggtitle(actual[counties[4], ]$name) +
    labs(col = NULL) +
    theme_bw()
  
  plt5 <- ggplot(data = county5, aes(x = years)) + 
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "pink", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
    geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
    xlab("year") + ylab("Relative Risk") + ggtitle(actual[counties[5], ]$name) +
    labs(col = NULL) +
    theme_bw()
  
  plt6 <- ggplot(data = county6, aes(x = years)) + 
    geom_ribbon(aes(x = years, ymin = lower_quant, ymax = upper_quant, col = "95% CI"), 
                fill = "pink", alpha = 0.6) +
    geom_line(aes(x = years, y = fitted_rate, col = "Fitted rate")) +
    geom_point(aes(x = years, y = true_rate, col = "True rate")) + 
    xlab("year") + ylab("Relative Risk") + ggtitle(actual[counties[6], ]$name) +
    labs(col = NULL) +
    theme_bw()
  
  ggarrange(plt1, plt2, plt3,
            plt4, plt5, plt6,
            ncol = 2, nrow = 3,
            common.legend = TRUE,
            legend = "right")
  
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






