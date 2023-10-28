
### Plotting functions 

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







