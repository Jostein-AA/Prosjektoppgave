#Clear environment
rm(list = ls())

#Load libraries
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

#Set working directory if not correct
if(getwd() != "C:/Users/joste/Documents/H2023/Code/Prosjektoppgave"){
  setwd("H2023/Code/Prosjektoppgave/")
}
source("utilities.R")

#Load data
ohio_df <- read.csv("ohio_df.csv")

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]

#Get number of counties (n) and number of years (T)
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points

###
#Plot mean and standard deviation of rate for each county

#Start by calculating mean and std. of rate for each county
ohio_map$mean_rate <- rep(0, n)  #Initialize to zero
ohio_map$sd_rate <- rep(0, n)  #Initialize to zero
for(county_name in ohio_df$name){ #Iterate over all the counties
  #Extract values for county 'county_name'
  temp <- ohio_df[ohio_df$name == county_name, ]
  
  #Calculate mean and standard deviation of rate * 100 000 for that county
  ohio_map$mean_rate[ohio_map$NAME == county_name] = mean(temp$rate * 1E5)
  ohio_map$sd_rate[ohio_map$NAME == county_name] = sd(temp$rate * 1E5)
}




scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
#Create heatmap of mean rate pr. 100 000 
mean_rate_plot <- ggplot(data = ohio_map) +  
  geom_sf(aes(fill = mean_rate), #Plots rate per 100 000
          alpha = 1,
          color="black") + ggtitle("Mean rate\n per 100 000") + 
  theme(plot.title = element_text(size = 15),
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

sd_rate_plot <- ggplot(data = ohio_map) +  
  geom_sf(aes(fill = sd_rate), #Plots cases per thousand
          alpha = 1,
          color="black") + ggtitle("Standard deviation of rate\n per 100 000") + 
  theme(plot.title = element_text(size = 15),
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

ggarrange(mean_rate_plot, sd_rate_plot,
          common.legend = FALSE)


###
#Heatmap plotted for each year

#Merge data with map for plotting over the areas
ohio_df_w_geometries <- merge(ohio_map, ohio_df,
                              by.x = c("NAME"), by.y = c("name"),
                              all = T, suffixes = T)

#Get rate per 100 000
ohio_df_w_geometries$rate <- ohio_df_w_geometries$rate * 1E5

min_rate <- min(ohio_df_w_geometries$rate); max_rate <- max(ohio_df_w_geometries$rate)
hardcoded_bins = seq(min_rate, max_rate, length.out = 8)
hardcoded_bins = round(hardcoded_bins, 0)

p_1 <- case_count_plot_1_year(ohio_df_w_geometries, 1968, hardcoded_bins = hardcoded_bins)
p_2 <- case_count_plot_1_year(ohio_df_w_geometries, 1969, hardcoded_bins = hardcoded_bins)
p_3 <- case_count_plot_1_year(ohio_df_w_geometries, 1970, hardcoded_bins = hardcoded_bins)
p_4 <- case_count_plot_1_year(ohio_df_w_geometries, 1971, hardcoded_bins = hardcoded_bins)
p_5 <- case_count_plot_1_year(ohio_df_w_geometries, 1972, hardcoded_bins = hardcoded_bins)
p_6 <- case_count_plot_1_year(ohio_df_w_geometries, 1973, hardcoded_bins = hardcoded_bins)
p_7 <- case_count_plot_1_year(ohio_df_w_geometries, 1974, hardcoded_bins = hardcoded_bins)
p_8 <- case_count_plot_1_year(ohio_df_w_geometries, 1975, hardcoded_bins = hardcoded_bins)
p_9 <- case_count_plot_1_year(ohio_df_w_geometries, 1976, hardcoded_bins = hardcoded_bins)
p_10 <- case_count_plot_1_year(ohio_df_w_geometries, 1977, hardcoded_bins = hardcoded_bins)
p_11 <- case_count_plot_1_year(ohio_df_w_geometries, 1978, hardcoded_bins = hardcoded_bins)
p_12 <- case_count_plot_1_year(ohio_df_w_geometries, 1979, hardcoded_bins = hardcoded_bins)
p_13 <- case_count_plot_1_year(ohio_df_w_geometries, 1980, hardcoded_bins = hardcoded_bins)
p_14 <- case_count_plot_1_year(ohio_df_w_geometries, 1981, hardcoded_bins = hardcoded_bins)
p_15 <- case_count_plot_1_year(ohio_df_w_geometries, 1982, hardcoded_bins = hardcoded_bins)
p_16 <- case_count_plot_1_year(ohio_df_w_geometries, 1983, hardcoded_bins = hardcoded_bins)
p_17 <- case_count_plot_1_year(ohio_df_w_geometries, 1984, hardcoded_bins = hardcoded_bins)
p_18 <- case_count_plot_1_year(ohio_df_w_geometries, 1985, hardcoded_bins = hardcoded_bins)
p_19 <- case_count_plot_1_year(ohio_df_w_geometries, 1986, hardcoded_bins = hardcoded_bins)
p_20 <- case_count_plot_1_year(ohio_df_w_geometries, 1987, hardcoded_bins = hardcoded_bins)
p_21 <- case_count_plot_1_year(ohio_df_w_geometries, 1988, hardcoded_bins = hardcoded_bins)


ggarrange(p_1, p_2, p_3, p_4, 
          p_5, p_6, p_7, p_8,
          p_9, p_10, p_11, p_12, 
          p_13, p_14, p_15, p_16,
          p_17, p_18, p_19, p_20,
          p_21,
          ncol = 5, nrow = 5, 
          common.legend = TRUE, legend = "left") + 
  theme(plot.margin = margin(0.1, 0, 0, 0, "cm")) 




###
#Depth measure
#Depth measure over time series

#Calculate rate per 100 000
ohio_df$rate <- ohio_df$rate * 1E5

#Transform data to functional data format (only for calculating depth measures)
grid = unique(ohio_df$year)
temp_ = ohio_df[ohio_df$name == unique(ohio_df$name)[1], ]
temp_data = temp_$rate
for(county in unique(ohio_df$name)[2:length(unique(ohio_df$name))]){
  temp_ = ohio_df[ohio_df$name == county, ]
  temp_data = rbind(temp_data, temp_$rate)
}

fdata <- fData(grid, temp_data)

#Calculate the median curve by way of modified band depth
median_curve = data.frame(
  median_curve = median_fData(fdata, type = "MBD")$values[1, ],
  year = grid)


#Plot the functional data (each county has its own curve)
#along with the depth median to see potential temporal trend in data
ohio_df$county <- ohio_df$name
theme_set(theme_bw())
ggplot(data = ohio_df, 
       aes(x = year, 
           y = rate, 
           group = county, 
           color = county),
       show.legend = FALSE) + 
  theme(legend.position="none") +
  geom_line() + 
  geom_line(data = median_curve, 
            aes(x = year, y = median_curve), 
            col = 'black', lwd = 3) +
  xlab("year") + ylab("rate per 100 000")














