
### Plotting functions 

#Plot heatmap for 1 year
case_count_plot_1_year <- function(sf_data, year, bins){
  scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
  scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
  #scale = heat.colors(8, rev= TRUE)
  ggplot(data = sf_data[sf_data$year == year, ]) + 
    geom_sf(aes(fill = rate), #Plots cases per thousand
            alpha = 1,
            color="black") + ggtitle(year) +
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
      breaks = bins,
      limits= c(0,100),
      guide = "colorscale")
}

#Plot heatmap for all years
plot_all_years <- function(sf_data, bins,
                           ncol = 7, nrow = 6,
                           widths =c(2, 2, 2, 2)
                           ){
  par(mar = c(0, 0, 0, 0))
  void_plot <- ggplot() + theme_void()
  #Hardcoded due to ggplot being a bitch
  par(mar = c(0, 0, 0, 0))
  plt1 <- case_count_plot_1_year(sf_data, 1, bins)
  plt2 <- case_count_plot_1_year(sf_data, 2, bins)
  plt3 <- case_count_plot_1_year(sf_data, 3, bins)
  plt4 <- case_count_plot_1_year(sf_data, 4, bins)
  plt5 <- case_count_plot_1_year(sf_data, 5, bins)
  plt6 <- case_count_plot_1_year(sf_data, 6, bins)
  plt7 <- case_count_plot_1_year(sf_data, 7, bins)
  plt8 <- case_count_plot_1_year(sf_data, 8, bins)
  plt9 <- case_count_plot_1_year(sf_data, 9, bins)
  plt10 <- case_count_plot_1_year(sf_data, 10, bins)
  plt11 <- case_count_plot_1_year(sf_data, 11, bins)
  plt12 <- case_count_plot_1_year(sf_data, 12, bins)
  plt13 <- case_count_plot_1_year(sf_data, 13, bins)
  plt14 <- case_count_plot_1_year(sf_data, 14, bins)
  plt15 <- case_count_plot_1_year(sf_data, 15, bins)
  plt16 <- case_count_plot_1_year(sf_data, 16, bins)
  plt17 <- case_count_plot_1_year(sf_data, 17, bins)
  plt18 <- case_count_plot_1_year(sf_data, 18, bins)
  plt19 <- case_count_plot_1_year(sf_data, 19, bins)
  plt20 <- case_count_plot_1_year(sf_data, 20, bins)
  plt21 <- case_count_plot_1_year(sf_data, 21, bins)
  
  
  ggarrange(plt1, void_plot, plt2, void_plot, plt3, void_plot, plt4,
            plt5, void_plot, plt6, void_plot, plt7, void_plot, plt8,
            plt9, void_plot, plt10, void_plot, plt11, void_plot, plt12,
            plt13, void_plot, plt14, void_plot, plt15, void_plot, plt16,
            plt17, void_plot, plt18, void_plot, plt19, void_plot, plt20,
            plt21,
               ncol = ncol, nrow = nrow, common.legend = TRUE, legend = "right",
            widths = widths
               )
}

### Model choice and summary functions







