#########################################################
## Creating the local region map/raster off the NC EEZ ##
#########################################################


library(devtools)
library(rgdal)
library(dplyr)
library(sp)
library(tools)
library(ncdf4)
library(raster)
library(rgeos)
library(plyr)

library(ggmap)
library(dplyr)
library(stringr)
library(readr)
library(pryr)
library(marmap)
library(lubridate)


# @knitr EARTHgrid

# Reading/creating the earth GRID 

earthGrid<-nc_open("./Environment/temp/earthgrid/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used

earthGrid<-ncvar_get(earthGrid,"thetao")

earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1 # 1 = earth

earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # basicly Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus

earthGrid <- reclassify(earthGrid, cbind(NA,2)) # 2 = ocean



# @knitr EEZgrid


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



# @knitr SubstrateCrop

# Crop the geomorphic substrate and coral


shps <- dir("./Environment/geomorph", "*.shp")
shps <- file_path_sans_ext(shps)

shps <- shps[-which(shps == "Abyssal_Classification" | shps == "Shelf_Classification")] # Removing the classification layers and keeping the simple versions of Abyss and Shelf

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


# Merging the features of the polygons to one polygon # Changed due to issues in the resulting polygons

croppedNcPoly <- croppedNcshp #lapply(croppedNcshp, gUnionCascaded)

croppedNcPoly <- croppedNcPoly[lapply(croppedNcPoly, length)!=0] # Removing potential null polygons


# Filtering the features through the occurrences and where they intersect (deleting features with no occurrences on them)

NCRecords <- read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE)

# keeping coords and species of interest
NCRecords <- NCRecords[,c(2:3, 7, 9, 11)]
NCRecords <- NCRecords[NCRecords$C_amblyrhynchos != 0 | NCRecords$C_melanopterus != 0 | NCRecords$T_obesus !=0,]
NCRecords$occurrence <- 1
rownames(NCRecords) <- NULL


## creating spatial points for our occurrences

RecPoints = SpatialPointsDataFrame(
  coords = cbind(NCRecords$longitude, NCRecords$latitude),
  data = NCRecords,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)

NCRecords <- NCRecords[which(raster::extract(eezNcGrid,RecPoints) == 1),]

RecPoints = SpatialPointsDataFrame(
  coords = cbind(NCRecords$longitude, NCRecords$latitude),
  data = NCRecords,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)


## Retrieving cells with occurrences in it

RecPoly <- raster(extent(eezNcGrid))
res(RecPoly) <- res(eezNcGrid)
projection(RecPoly) <- proj4string(eezNcGrid)
origin(RecPoly) <- origin(eezNcGrid)
RecPoly <- rasterize(RecPoints,field = "occurrence", RecPoly)

RecPoly <- as(RecPoly, "SpatialPolygons") # Points are now one with cell under a polygon form



# retrieving features that intersect with the occurrence cell

crossingFeatures <- lapply(croppedNcPoly, function(x, Rec = RecPoly){ 
  
  true = which(gIntersects(x,RecPoly, byid = TRUE) == TRUE)
  
  true
  
})



# Intersecting the feature's polygons with the occurrences


crpdNcPolyFiltrd <- croppedNcPoly[lapply(crossingFeatures, length)!=0]


interCrss <- lapply(crpdNcPolyFiltrd,function(x, e = eezNcPolyGrid){
  
  inter = gIntersection(x, eezNcPolyGrid, byid = TRUE)
  
  inter
})


# @knitr ExtractCov 

# Creating the final data.frame with cells numbers and values associated

gridCells <- as.data.frame(eezNcGrid)




# extract percentages of each polygons in each cells

getCellsDf <- function(interCrss) {


  featCellsCoords <- lapply(interCrss, coordinates)

  featCells <- lapply(featCellsCoords, extract, x = eezNcGrid,cellnumber = TRUE)

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

# substrateLs <- lapply(featCov, getSingleGridData) # If getSingleGridData returned abyss
# 
# substrate <- as.data.frame(substrateLs) # If getSingleGridData returned abyss
# 
# colnames(substrate) <- names(substrateLs) # If getSingleGridData returned abyss


  



