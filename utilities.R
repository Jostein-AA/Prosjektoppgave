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


#Plot posterior intercept
plot_intercept <- function(fitted_model){
  #Function that plots the posterior distribution of the intercept
  par(mfrow = c(1, 1))
  plot(fitted_model$marginals.fixed$`(Intercept)`[, 1],
       fitted_model$marginals.fixed$`(Intercept)`[, 2],
       xlab = "", ylab = "", type = "l", lwd = 2.5, 
       main = "Posterior density of Intercept")
  
}

#Plot posterior precision's
plot_precisions_random_effects <- function(fitted_model){
  par(mfrow = c(2, 2))
  #Plot precision of iid temporal effect
  plot(fitted_model$marginals.hyperpar$`Precision for year (iid component)`[, 1],
       fitted_model$marginals.hyperpar$`Precision for year (iid component)`[, 2],
       type = "l", xlab = "", ylab = "", main = "Precision of iid temporal effect")
  #Plot precision of structured temporal effect
  plot(fitted_model$marginals.hyperpar$`Precision for year (spatial component)`[, 1],
       fitted_model$marginals.hyperpar$`Precision for year (spatial component)`[, 2],
       type = "l", xlab = "", ylab = "", main = "Precision of structured temporal effect")
  #Plot precision of iid spatial effect
  plot(fitted_model$marginals.hyperpar$`Precision for county (iid component)`[, 1],
       fitted_model$marginals.hyperpar$`Precision for county (iid component)`[, 2],
       type = "l", xlab = "", ylab = "", main = "Precision of iid spatial effect")
  #Plot precision of structured spatial effect
  plot(fitted_model$marginals.hyperpar$`Precision for county (spatial component)`[, 1],
       fitted_model$marginals.hyperpar$`Precision for county (spatial component)`[, 2],
       type = "l", xlab = "", ylab = "", main = "Precision of structured spatial effect")
}


#Plot the temporal effect
plot_temporal_effect <- function(fitted_model, T){
  par(mfrow = c(1, 1))
  matplot(fitted_model$summary.random$year[(T + 1):(2 * T), 4:6],
          lty=c(2,1,2), type="l", col=1,
          xlab = "year", ylab = "Temporal random effect")
}

#Plot the spatial effect as a heatmap
plot_spatial_effect <- function(map, fitted_model, n){
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale_1 = scale_col[c(3,10,13,17,21,24,27,30)] #Select color scale to be more red
  
  spatial_structured_effect_mean <- fitted_model$summary.random$county$mean[(n+1):(2*n)]
  spatial_structured_effect_sd <- fitted_model$summary.random$county$sd[(n+1):(2*n)]
  temp_ohio_map <- map[ ,c("geometry", "NAME")]
  temp_ohio_map$mean <- spatial_structured_effect_mean
  temp_ohio_map$sd <- spatial_structured_effect_sd
  
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
            color="black") + ggtitle("Std deviation of spatial effect for each county") +
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



#Plot precision of interaction
plot_prec_interactions <- function(fitted_model, title){
  par(mfrow = c(1, 1))
  plot(fitted_model$marginals.hyperpar$`Precision for space.time`[, 1],
       fitted_model$marginals.hyperpar$`Precision for space.time`[, 2],
       type = "l", xlab = "", ylab = "", lwd = 2.5,
       main = title)
}

#Plot the interactions
plot_interaction <- function(fitted_model){
  matplot(fitted_model$summary.random$space.time[ , 4:6],
          lty=c(2,1,2), type="l", col=1,
          xlab = "Datum ID", ylab = "Interaction effect")
}


#Plot fitted values against actual values (all together)
plot_fitted_vs_actual_together <- function(actual_data, fitted_model){
  par(mfrow = c(1, 1))
  plot(actual_data$rate,
       xlab = "Datum ID", ylab = "Fitted value") + 
    lines(fitted_model$summary.fitted.values[4],
          lty = 1, type = "l", 
          col = 1)
}

#Plot a time series for every county showing fitted values vs actual values
every_county_time_series <- function(fitted_model, n, T){
  #model_fitted: a model that is fitted, from which fitted values can be extracted
  #Plots a time
  par(mfrow=c(3,3))
  for(i in 1:n){
    plot(ohio_df$rate[seq(i,n*T,by=n)], 
         ylab = "rate",
         xlab = "year",
         main = ohio_df$county_name[i])
    matplot(fitted_model$summary.fitted.values[seq(i,n*T,by=n),3:5],
            col=1, lty=c(2,1,2), type="l", add=T)
  }
}







#Plot heatmap for 1 year
case_count_plot_1_year <- function(sf_data, year, hardcoded_bins = NULL){
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
  #scale = heat.colors(8, rev= TRUE)
  
  if(is.null(hardcoded_bins)){
    p <- ggplot(data = sf_data[sf_data$year == year, ]) + 
      geom_sf(aes(fill = rate), #Plots death rate
              alpha = 1,
              color="black") + 
      ggtitle(year) +
      theme(plot.title = element_text(size = 5),
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
  } else{
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
  }
  return(p)
}


#Plot heatmap for all years
plot_all_years <- function(sf_data){
  #Hardcoded due to ggplot 
  plots = c()
  plots[1] <- case_count_plot_1_year(sf_data, 1)
  plots[2] <- case_count_plot_1_year(sf_data, 2)
  plots[3] <- case_count_plot_1_year(sf_data, 3)
  plots[4] <- case_count_plot_1_year(sf_data, 4)
  plots[5] <- case_count_plot_1_year(sf_data, 5)
  plots[6] <- case_count_plot_1_year(sf_data, 6)
  plots[7] <- case_count_plot_1_year(sf_data, 7)
  plots[8] <- case_count_plot_1_year(sf_data, 8)
  plots[9] <- case_count_plot_1_year(sf_data, 9)
  plots[10] <- case_count_plot_1_year(sf_data, 10)
  plots[11] <- case_count_plot_1_year(sf_data, 11)
  plots[12] <- case_count_plot_1_year(sf_data, 12)
  plots[13] <- case_count_plot_1_year(sf_data, 13)
  plots[14] <- case_count_plot_1_year(sf_data, 14)
  plots[15] <- case_count_plot_1_year(sf_data, 15)
  plots[16] <- case_count_plot_1_year(sf_data, 16)
  plots[17] <- case_count_plot_1_year(sf_data, 17)
  plots[18] <- case_count_plot_1_year(sf_data, 18)
  plots[19] <- case_count_plot_1_year(sf_data, 19)
  plots[20] <- case_count_plot_1_year(sf_data, 20)
  plots[21] <- case_count_plot_1_year(sf_data, 21)
  
  
  #ggarrange(plt1, void_plot, plt2, void_plot, plt3, void_plot, plt4,
  #          plt5, void_plot, plt6, void_plot, plt7, void_plot, plt8,
  #          plt9, void_plot, plt10, void_plot, plt11, void_plot, plt12,
  #          plt13, void_plot, plt14, void_plot, plt15, void_plot, plt16,
  #          plt17, void_plot, plt18, void_plot, plt19, void_plot, plt20,
  #          plt21,
  #             ncol = ncol, nrow = nrow, common.legend = TRUE, legend = "right")
  return(plots)
}

### Model choice and summary functions







