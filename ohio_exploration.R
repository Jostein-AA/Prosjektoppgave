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

#Load data
ohio_df <- read.csv("ohio_df.csv")

#read the shapefile of ohio
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]


#Start by finding median/mean/sd rate for each county
ohio_map$median_rate <- rep(0, n)  #Initialize to zero
ohio_map$mean_rate <- rep(0, n)  #Initialize to zero
ohio_map$sd_rate <- rep(0, n)  #Initialize to zero
for(county_name in ohio_df$county_name){ #Iterate over all the counties
  #Extract values for county 'county_name'
  temp <- ohio_df[ohio_df$county_name == county_name, ]
  
  #Calculate median/mean/sd rate for that county
  ohio_map$median_rate[ohio_map$NAME == county_name] = median(temp$rate)
  ohio_map$mean_rate[ohio_map$NAME == county_name] = mean(temp$rate)
  ohio_map$sd_rate[ohio_map$NAME == county_name] = sd(temp$rate)
}

scale_col = heat.colors(30, rev=TRUE) #Divide color gradient into 30 
scale = scale_col[c(3,10,13,18,21,24,27,30)] #Select color scale to be more red
#scale = heat.colors(8, rev= TRUE)
mean_rate_plot <- ggplot(data = ohio_map) +  
  geom_sf(aes(fill = mean_rate), #Plots cases per thousand
          alpha = 1,
          color="black") + ggtitle("Mean rate each county") + 
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
    palette = function(x) c(scale),
    labels = function(x){x},
    guide = "colorscale")

sd_rate_plot <- ggplot(data = ohio_map) +  
  geom_sf(aes(fill = sd_rate), #Plots cases per thousand
          alpha = 1,
          color="black") + ggtitle("Standard deviation of rate\n for each county") + 
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
    palette = function(x) c(scale),
    labels = function(x){x},
    guide = "colorscale")

ggarrange(mean_rate_plot, sd_rate_plot,
          common.legend = FALSE)





















