library(spData)
library(sf)
library(spdep)
library(ggplot2)

# read the shapefile, defining the country and subregional borders
map <- st_read(system.file("shapes/columbus.shp",
                           package = "spData"), quiet = TRUE)


# extract the neighbourhood structure from the graph
nb <- spdep::poly2nb(map, queen = FALSE)
head(nb)

# plot the country indication the first-order neighbours
plot(st_geometry(map), border = "lightgrey")
plot.nb(nb, st_geometry(map), add = TRUE)

# needed for inla
file4inla <- nb2INLA("myadj.txt", nb)
# or
matrix4inla <- nb2mat(nb, style="B")
mydiag = rowSums(matrix4inla)
matrix4inla <- -matrix4inla
diag(matrix4inla) <- mydiag

# plot neighbours for region 20. 
id <- 20 # area id
map$neighbors <- "other"
map$neighbors[id] <- "area"
map$neighbors[nb[[id]]] <- "neighbors"
ggplot(map) + geom_sf(aes(fill = neighbors)) + theme_bw() +
  scale_fill_manual(values = c("gray30", "gray", "white"))