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

#Load other files
source("./functions_plotting_etc.R")


#Load data: contains deaths, county (not named), year, and expected deaths
#raw_df <- read.table("http://faculty.washington.edu/jonno/SISMIDmaterial/ohiodata2_spacetime.txt",header=T)

#Load second dataset, contains: county, gender, race, year, deaths, number at risk, county name
raw_df_2 <- read.csv("dataohiocomplete.csv")

# read the shapefile of Ohio and sort rows alphabeticaly on county names
ohio_map  <- read_sf('./fe_2007_39_county')[ ,c("geometry", "NAME")]
ohio_map <- ohio_map[order(ohio_map$NAME), ]


#Transform year from 1968,... to 1,...
raw_df_2$year <- raw_df_2$year - min(raw_df_2$year) + 1

#Rename
colnames(raw_df_2) <- c("county", "gender", "race",
                        "year", "observed", "pop_at_risk", "county_name") 

#Aggregate observed deaths over gender and race
observed_aggregated <- aggregate(x = raw_df_2$observed, by = list(raw_df_2$county,
                                                                  raw_df_2$county_name,
                                                                  raw_df_2$year), FUN = sum)
#Aggregate pop_at_risk over gender and race
pop_at_risk_aggregated <- aggregate(x = raw_df_2$pop_at_risk, by = list(raw_df_2$county,
                                                                        raw_df_2$county_name,
                                                                        raw_df_2$year), FUN = sum)


#Join observed and expected aggregate on Group.3 (year), group.2 (county_name), group.1 (county)
ohio_df <- merge(observed_aggregated, pop_at_risk_aggregated,
                 by = c("Group.3", "Group.2", "Group.1"),
                 all = T, suffixes = T)

#Sort years in increasing order (if not, years are not sorted)
ohio_df = ohio_df[order(ohio_df$Group.3, decreasing = F), ]

#Rename to understandable column names.
colnames(ohio_df) <- c("year", "county_name", "county", "deaths", "pop_at_risk")

#Check that ohio_df$observed correspond with raw_df$Y (should be the same)
#temp <- ohio_df$deaths - raw_df$Y
#if(length(unique(temp)) == 1 & temp[1] == 0){
#  print("They correspond")
  
  #If they correspond, give ohio_df the expected number of deaths for each county each year
#  ohio_df$expected <- raw_df$E
#} else{
#  print("They do NOT correspond")
#}

#Create copies of time for unstructured/structured random effects
ohio_df$time.unstructured <- ohio_df$time.structured <- ohio_df$year

#Get death rate instead (NB: Only for plotting)
ohio_df$rate <- ohio_df$deaths/ohio_df$pop_at_risk

#Add constant to ohio_df
#ohio_df$dummy <- rep(1, length(ohio_df$year))


#Find dimensions
n <- length(unique(ohio_df$county))   # Number of areas
T <- length(unique(ohio_df$year))     # Number of time points

#structure matrices for random walk 1 and 2
struct_RW1 <- INLA:::inla.rw(n = T, order = 1, 
                             scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                             sparse = TRUE)

struct_RW2 <- INLA:::inla.rw(n = T, order = 2, 
                             scale.model = FALSE, # set scale.model  = F because we'll scale in the formula
                             sparse = TRUE)

#Create ICAR structure matrix
nb <- spdep::poly2nb(ohio_map, queen = FALSE)
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag
ICAR_structure <- matrix4inla

# quick check: If you coded your ICAR precision matrix correctly, both rows and
# columns will sum to zero
print("rows sum-to-zero:")
print((length(unique(apply(ICAR_structure, 1, sum))) == 1 &
         apply(ICAR_structure, 1, sum)[1] == 0))
print("columns sum-to-zero:")
print((length(unique(apply(ICAR_structure, 2, sum))) == 1 &
         apply(ICAR_structure, 2, sum)[1] == 0))

#Make the ICAR_structure matrix sparse.
ICAR_structure <- Matrix(ICAR_structure, sparse = TRUE) #Make it sparse


#Add space_time unstructured to ohio_df
ohio_df$space_time_unstructured <- 1:(n *T)

#Save data and precision matrices
write.csv(ohio_df, file = "ohio_df.csv", row.names = FALSE)

writeMM(struct_RW1, file = 'struct_RW1.txt')
writeMM(struct_RW2, file = 'struct_RW2.txt')
writeMM(ICAR_structure, file = 'ICAR_struct.txt')

#How to read the precision matrices
struct_RW1.test = readMM(file='struct_RW1.txt')
struct_RW2.test = readMM(file = 'struct_RW2.txt')
ICAR_structure.test = readMM(file = 'ICAR_struct.txt')

ohio_df.test <- read.csv("ohio_df.csv")


#Sanity check on ICAR
#-------------------------------------------------------------------------------
#Check that adjacencies are actually adjacencies
plot(st_geometry(ohio_map), border = "lightgrey")
plot.nb(nb, st_geometry(ohio_map), add = TRUE)

#How to check that the correct counties border the correct counties???
#If it does not print, it means the Names of the counties line up
for(i in 1:length(ohio_map$NAME)){
  if(ohio_map$NAME[i] != ohio_df$county_name[i]){
    print("county names do not correspond")
  }
}


#Check that county id's in ohio_df correspond to ohio_map index
#If no print of Uh oh, its ok
for(i in 1:length(ohio_map$NAME)){
  if(i != ohio_df$county[i]){
    print("County id's and county_map index do not correspond")
  }
}

#Have to check that ICAR_structure matches neighborhood structure
equal = 0
for(i in 1:n){
  temp = which(ICAR_structure[i, ] != 0)
  temp = temp[which(temp != i)]
  
  #The ICAR is based on the 
  if(identical(temp, nb[[i]])){
    equal = equal + 1
  }
}
print("Number of equals (should be 88)")
print(equal)

#It should all match up


map = ohio_map

# Check that for area (id) that neighbors are adjacent 
id <- 53 # area id
map$neighbors <- "other"
map$neighbors[id] <- "area"
map$neighbors[nb[[id]]] <- "neighbors"
ggplot(map) + geom_sf(aes(fill = neighbors)) + theme_bw() +
  scale_fill_manual(values = c("gray30", "gray", "white"))
#-------------------------------------------------------------------------------



















