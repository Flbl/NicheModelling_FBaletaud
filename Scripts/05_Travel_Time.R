#################################################################
##                                                             ##
##          Getting and formatting travel time data            ##         
##                                                             ##
#################################################################

# Transforming points coordinates to Behrman


eezNcGrid <- raster("./data/interdata/eezNcGrid.tif")

noumea <- as.matrix(data.frame(x = 166.4392, y = -22.2558)) # Arbitrary coordinates of noumea that fall in cell where noumea is
eezNcGrid[extract(eezNcGrid, noumea, cellnumbers = TRUE)[1]] <- 1


ncGridBehr <- projectRaster(eezNcGrid, 
                            crs = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"), 
                            method = "ngb",
                            res = 4000)

plot(ncGridBehr)

plot(eezNcGrid)

writeRaster(ncGridBehr, "./data/interdata/eezNcGridBehrman.tif", overwrite = TRUE) 

noumeaBehr <- SpatialPoints(coords = noumea,
                            proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
noumeaBehr <- spTransform(noumeaBehr, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))
coordinates(noumeaBehr) # Coordinates of noumea in behrman projection to enter in the python function args


# Calling python function to create the raster of distance 








# reading the generated raster


costNoumea <- raster("./data/interdata/costNoumea_0.img", crs = CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs") )
costNoumea
plot(costNoumea)

costNoumeaLatLong <- projectRaster(from = costNoumea, to = eezNcGrid, 
                                   method = "ngb")

# plot(costNoumeaLatLong)

writeRaster(costNoumeaLatLong, filename = "./data/predictdata/Travel_Dist.tif", overwrite = TRUE)







