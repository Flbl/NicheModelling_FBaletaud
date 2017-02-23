#################################################################
##                                                             ##
##  Getting and formatting the geomorphic features variables   ##         
##                        (substrate)                          ##
#################################################################

library(devtools)
library(rgdal)
library(plyr)
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



# @knitr EEZgrid

eezNcGrid <- raster("./data/interdata/eezNcGrid.tif")
eezNcPoly <- readOGR(dsn = "./data/interdata", layer = "eezNcPoly")


# Converting to polygons 

eezNcPolyGrid <- as(eezNcGrid, "SpatialPolygons")

eezNcPolyGrid$cellID <- c(1:length(eezNcPolyGrid))



# @knitr SubstrateCrop

# Crop the geomorphic substrate and coral


shps <- dir("./data/rawdata/Environment/geomorph", "*.shp")
shps <- file_path_sans_ext(shps)

shps <- shps[which(shps == "Shelf" | shps == "Slope" | shps == "Abyss"| shps == "Canyons"| shps == "Guyots"| shps == "Seamounts"| shps == "Ridges"| shps == "Escarpments")]

# shps <- shps[-which(shps == "Abyssal_Classification" | shps == "Shelf_Classification")] # Removing the classification layers and keeping the simple versions of Abyss and Shelf

# Function to crop 

cropNc <- function(shp, cropNc = eezNcPoly, PolyGrid = eezNcPolyGrid ){
  
  shp <- assign(shp, readOGR("./data/rawdata/Environment/geomorph",layer = shp))
  
  shp <- crop(shp, extent(cropNc))
  
  shp
}


croppedNcshp <- lapply(shps, cropNc)

names(croppedNcshp) <- shps

# croppedNcshp <- croppedNcshp[lapply(croppedNcshp, length)!=0]


#Adding the coral

## Crop the Coral 

CoralShp = readOGR("./data/rawdata/Environment/WCMC008_CoralReef2010_v1_3/01_Data", layer = "14_001_WCMC008_CoralReef2010_v1_3")

CoralNc <- crop(CoralShp, extent(eezNcPoly))


croppedNcshp$Coral <- CoralNc


# Merging the features of the polygons to one polygon # Changed due to issues in the resulting polygons

croppedNcPoly <- croppedNcshp # before lapply(croppedNcshp, gUnionCascaded)

croppedNcPoly <- croppedNcPoly[lapply(croppedNcPoly, length)!=0] # Removing potential null polygons

# Intersecting the features with each "cellpolygons" of eezNcPolyGrid to get a polygon of the features for each cell it intersects with

interCrss <- lapply(croppedNcPoly,function(x, e = eezNcPolyGrid){
  
  inter = gIntersection(x, e, byid = TRUE)
  
  inter
})


# @knitr ExtractCov 

# Creating the final data.frame with cells numbers and values associated

gridCells <- as.data.frame(eezNcGrid)
colnames(gridCells) <- "layer"



# extract percentages of each polygons in each cells

getCellsDf <- function(interCrss) {


  featCellsCoords <- lapply(interCrss, coordinates)

  featCells <- lapply(X = featCellsCoords, FUN = extract, x = eezNcGrid, cellnumber = TRUE)

  featCells <- lapply(featCells, as.data.frame)
  
  featCells

}# eo getCellsDf




featCells <- getCellsDf(interCrss)



behrSubst <- lapply(interCrss, spTransform, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")) 



getSurfPoly <- function(behrSubst){
  
  featCov <- lapply(behrSubst, gArea, byid = TRUE)
  
  featCov <- lapply(featCov, function(x){
    
    x <- x / (1000*1000)
  
    x
})
 
  featCov 
  
}# eo getSurfPoly



featCov <- getSurfPoly(behrSubst)

featCov <- mapply(cbind, featCells, area = featCov, SIMPLIFY = FALSE)




getSingleGridData <- function(x, gridd = gridCells){
  
  substrate <- gridd
  
  substrateNoNa <- data.frame(layer = gridd$layer[!is.na(gridd$layer)])
  
  rownames(substrateNoNa) <- na.omit(rownames(as.data.frame(substrate))[substrate[["layer"]] == 1])
  
  substrateNoNa$RN <- as.numeric(rownames(substrateNoNa))
  
  abyss1 <- data.frame(area = x[,3], RN = x[,1])
  
  abyss <- aggregate(abyss1, by = list(cells = abyss1$RN), sum) # some polygons are in the same cell so we add their surfaces
  
  rownames(abyss) <- abyss$cells
  
  zeroNames <- substrateNoNa$RN[is.na(match(substrateNoNa$RN, abyss$cells))]
  
  zeroDat <- data.frame(area = numeric(dim(substrateNoNa)[1] - dim(abyss)[1]), row.names = zeroNames)
  
  abyss <- data.frame(area = abyss$area, row.names = rownames(abyss))
  
  abyss <- rbind(abyss, zeroDat)
  
  abyss <- data.frame(area = abyss[ order(as.numeric(row.names(abyss))), ])
  
  rownames(abyss) <- rownames(substrateNoNa)
  
  # abyss # To directly get the dataframe of 1 column with the area of cells "in ocean" detailed in rownames (length = 17105)
  
  abyssNa <- rep (NA, dim(substrateNoNa)[1])

  abyssNa[!is.na(gridd)] <- abyss[,"area"]

  abyssNa # To get a complete list of the cells from eezNcGrid with the value associated (length = 29592)
  
}#eo getSingleGridData




substrate <- lapply(featCov, getSingleGridData)
substrate <- as.data.frame(substrate)
substrate <- cbind(coordinates(eezNcGrid), substrate)
substrate <- na.omit(substrate)


# @knitr WriteToRast

writeVarToRast <- function(s, map = eezNcGrid){
  
  sNoCoords <-  subset(s, select = -c(x,y))
  
  substrateSpatial <-  SpatialPointsDataFrame(
    coords = cbind(s$x, s$y),
    data = sNoCoords,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  featRast <- raster(extent(map))
  res(featRast) <- res(map)
  projection(featRast) <- proj4string(map)
  origin(featRast) <- origin(map)
  
  lapply(names(sNoCoords), function(x){
    
    rast <- featRast
    
    rast <- rasterize(substrateSpatial ,field = x, featRast)
    
    writeRaster(rast, filename = paste0("./data/calibdata/", x, ".tif"), overwrite = TRUE)
    
    
    
    
  })
    
    
  
}# eo writeVarToRast
    
writeVarToRast(substrate)

# plot(raster("./data/calibdata/Seamounts.tif"))










