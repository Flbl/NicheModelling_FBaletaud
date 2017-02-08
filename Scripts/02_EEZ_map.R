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
library(tools)
library(ncdf4)
library(raster)
library(rgeos)

library(ggmap)
library(dplyr)
library(stringr)
library(readr)
library(pryr)
library(marmap)
library(lubridate)
library(broom)
library(mregions)



# Reading/creating the earth GRID 

earthGrid<-nc_open("./Environment/temp/earthgrid/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used

earthGrid<-ncvar_get(earthGrid,"thetao")

earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1 # 1 = earth

earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # basicly Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus

earthGrid <- reclassify(earthGrid, cbind(NA,2)) # 2 = ocean






# read the eez shapefile and crop the NC eez

eez <- readOGR(dsn = "./Environment/World_EEZ_v9_20161021", layer = "eez")

# The 2 polygons of interest representing the NC EEZ
eez[eez$PolygonID == 31,]
eez[eez$PolygonID == 254,]

eezNcPoly <- eez[which(eez$PolygonID == 254 | eez$PolygonID == 31),]




# intersect earthGrid with only the New Caledonian polygons

eezNcGrid <- mask(earthGrid, eezNcPoly)
eezNcGrid <- crop(eezNcGrid, extent(eezNcPoly))
eezNcGrid <- reclassify(eezNcGrid, cbind(1,NA)) # NA = earth/outlimits
eezNcGrid <- reclassify(eezNcGrid, cbind(2,1)) # 1 = Ocean


# Converting to polygons 

eezNcPolyGrid <- as(eezNcGrid, "SpatialPolygons")

eezNcPolyGrid$cellID <- c(1:length(eezNcPolyGrid))





# Crop the geomorphic substrate and coral

## @knitr seafloor


shps <- dir("./Environment/geomorph", "*.shp")
shps <- file_path_sans_ext(shps)

shps <- shps[-c(2,17)] # Removing the classification layers and keeping the simple versions of Abyss and Shelf

# Function to crop 

cropNc <- function(shp, cropNc = eezNcPoly, PolyGrid = eezNcPolyGrid ){
  
  shp <- assign(shp, readOGR("./Environment/geomorph",layer = shp))
  
  shp <- crop(shp, extent(cropNc))
  
  shp
}


croppedNcshp <- lapply(shps, cropNc)

names(croppedNcshp) <- shps

croppedNcshp <- croppedNcshp[lapply(croppedNcshp, length)!=0]


#Adding the coral

## Crop the Coral 

CoralShp = readOGR("./Environment/WCMC008_CoralReef2010_v1_3/01_Data", layer = "14_001_WCMC008_CoralReef2010_v1_3")

CoralNc <- crop(CoralShp, extent(eezNcPoly))


croppedNcshp$Coral <- CoralNc


# Merging the features of the polygons to one polygon

cropNcPoly <- function(shp, polyGrid = eezNcPolyGrid){
  
  shp <- gUnionCascaded(shp)
  
  shp <- gIntersection(eezNcPolyGrid, shp, byid = TRUE)
  
  shp
  
}

croppedNcPoly <- lapply(croppedNcshp, cropNcPoly)

croppedNcPoly <- croppedNcPoly[lapply(croppedNcPoly, length)!=0]





# extract percentages of each polygons in each cells

#Created polygons have an area saved
#Need to find a way to associate the area of each polygon corresponding to a cell with the good cell ID from eezNcGridPoly created 



