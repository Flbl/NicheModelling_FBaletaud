#################################################################
##                                                             ##
##            Getting and formatting project Grids             ##         
##                                                             ##
#################################################################

library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(ncdf4)


# @knitr EARTHgrid

# Reading/creating the earth GRID 

earthGrid<-nc_open("./data/rawdata/Environment/temp/earthgrid/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used


earthGrid<-ncvar_get(earthGrid,"thetao")


earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1
earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus
earthGrid <- reclassify(earthGrid, cbind(NA,2)) # 2 = ocean

## Saving
writeRaster(earthGrid, "./data/interdata/earthGrid.tif", overwrite = TRUE)



# @knitr EEZgrid

# Creating the local region map/raster off the NC EEZ

## read the eez shapefile and crop the NC eez

eez <- readOGR(dsn = "./data/rawdata/Environment/World_EEZ_v9_20161021", layer = "eez")



## The 2 polygons of interest representing the NC EEZ

eez[eez$PolygonID == 31,]
eez[eez$PolygonID == 254,]

eezNcPoly <- eez[which(eez$PolygonID == 254 | eez$PolygonID == 31),]

## saving polygons representing the eez of New Caledonia
writeOGR(eezNcPoly, dsn = "./data/interdata", layer = "eezNcPoly", driver = "ESRI Shapefile")


## intersect earthGrid with only the New Caledonian polygons

eezNcGrid <- mask(earthGrid, eezNcPoly)
eezNcGrid <- crop(eezNcGrid, extent(eezNcPoly))
eezNcGrid <- reclassify(eezNcGrid, cbind(1,NA)) # NA = earth/outlimits
eezNcGrid <- reclassify(eezNcGrid, cbind(2,1)) # 1 = Ocean

## Saving
writeRaster(eezNcGrid, "./data/interdata/eezNcGrid.tif", overwrite = TRUE)




