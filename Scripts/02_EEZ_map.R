#########################################################
## Creating the local region map/raster off the NC EEZ ##
#########################################################

library(rgbif)
library(taxize)
library(taxizesoap)
library(devtools)
library(robis)
library(leaflet)
library(tibble)
library(rgdal)
library(dplyr)
library(sp)

library(ggmap)
library(dplyr)
library(ncdf4)
library(raster)
library(stringr)
library(readr)
library(pryr)
library(marmap)
library(lubridate)
library(broom)
library(mregions)
library(rgeos)


# Reading/creating the earth GRID a transformer en fonction

earthGrid<-nc_open("./Environment/temp/earthgrid/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used


earthGrid<-ncvar_get(earthGrid,"thetao")



earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1



earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus

earthGrid <- reclassify(earthGrid, cbind(NA,2))


earthRegion <- readOGR(dsn = "./Environment/GSHHS_region", layer = "GSHHS_f_L1")




# read the eez shapefile and crop the NC eez


eez <- readOGR(dsn = "./Environment/World_EEZ_v9_20161021", layer = "eez")

# The 2 polygons of interest representing the NC EEZ
eez[eez$PolygonID == 31,]
eez[eez$PolygonID == 254,]

eezNc1 <- eezNc[which(eez$PolygonID == 254 | eez$PolygonID == 31),]
plot(earthGrid)
plot(eezNc1, add = TRUE)

# gintersect with only the New Caledonian polygons
eezNc2 <- mask(earthGrid, eezNc1)
eezNc2 <- crop(eezNc2, extent(eezNc1))


