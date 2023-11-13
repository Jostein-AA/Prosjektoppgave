#Load libraries
library(INLA)
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
library(shapefiles)

#Set working directory
if(getwd() != "C:/Users/joste/Documents/H2023/Code/Prosjektoppgave"){
  setwd("H2023/Code/Prosjektoppgave/")
}

#Load raw dataset, contains: county, gender, race, year, deaths, number at risk, county name
raw_df <- read.csv("dataohiocomplete.csv")

# read the shapefile of Ohio and sort rows alphabeticaly on county names
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]


#Aggregate observed deaths over gender and race
deaths_aggregated <- aggregate(x = raw_df$y, by = list(raw_df$county,
                                                         raw_df$NAME,
                                                         raw_df$year), 
                                 FUN = sum)
#Aggregate pop_at_risk over gender and race
pop_at_risk_aggregated <- aggregate(x = raw_df$n, by = list(raw_df$county,
                                                            raw_df$NAME,
                                                            raw_df$year), 
                                    FUN = sum)


#Join observed and expected aggregate on Group.3 (year), group.2 (NAME), group.1 (county)
ohio_df <- merge(deaths_aggregated, pop_at_risk_aggregated,
                 by = c("Group.3", "Group.2", "Group.1"),
                 all = T, suffixes = T)

#Make sure years are sorted in increasing order
ohio_df = ohio_df[order(ohio_df$Group.3, decreasing = F), ]

#Rename to understandable column names.
colnames(ohio_df) <- c("year", "county_name", "county", "deaths", "pop_at_risk")


#Create copies of time for unstructured/structured random effects
#ohio_df$time.unstructured <- ohio_df$time.structured <- ohio_df$year

#Get death rate instead (NB: Only for plotting)
ohio_df$rate <- ohio_df$deaths/ohio_df$pop_at_risk


#Find dimensions
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points


#Add space.time for interactions to ohio_df
ohio_df$space.time <- 1:(n *T)

#Save data and precision matrices
write.csv(ohio_df, 
          file = "ohio_df.csv", 
          row.names = FALSE)


#Check that it is correct
ohio_df.test <- read.csv("ohio_df.csv")




















