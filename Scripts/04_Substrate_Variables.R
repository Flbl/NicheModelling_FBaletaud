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
    
    writeRaster(rast, filename = paste0("./data/predictdata/", x, ".tif"), overwrite = TRUE)
    
    
    
    
  })
    
    
  
}# eo writeVarToRast
    
writeVarToRast(substrate)

# plot(raster("./data/calibdata/Seamounts.tif"))


# @knitr ExctractCalibData



############ Extracting calib data ################## 

rawSubs <- substrate
rawSubs$cellNumber <- rownames(rawSubs)
head(rawSubs)


NCRecords <- read.table("./data/rawdata/Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE) 

species <- names(cleandOccs)

spDatasets <- lapply(species, function(my_sp, rec = NCRecords){

  
  matchSp = paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, "_"))[2], sep = "", collapse = "")
  matchCol = grep(matchSp,names(rec))
  names(rec)[matchCol] <- my_sp
  
  NC_occs <- rec[,c(3,2, matchCol)]
  NC_occs <- NC_occs[NC_occs[,3] >= "1",]
  colnames(NC_occs) <- c("decimalLongitude","decimalLatitude","individualCount")
  NC_occs$occurrence = 1
  NC_occs$individualCount <- NULL
  rownames(NC_occs) <- NULL # Reassign the rownames of NC_occs
  
  NC_occs$species <- my_sp
  
  NC_occs <- NC_occs[,c(4,1,2,3)]
  
  NC_occs$cellNumber <- cellFromXY(eezNcGrid, as.matrix(NC_occs[,c(2,3)]))
  
  ## Checking for points on earth
  ## creating spatial points for our occurrences
  
  occsSpatial = SpatialPointsDataFrame(
    coords = as.matrix(NC_occs[,c(2,3)]),
    data = NC_occs,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  )
  
  NC_occs <- NC_occs[-which(is.na(raster::extract(eezNcGrid,occsSpatial))),]
  
  
  NC_occs
  
})

names(spDatasets) <- species

# spDatasets

## Generating absences

generateAbs <- function(Occs){
  
  # Creating pseudo absences using IUCN range map outer limits
  print("Getting IUCN range map data")
  
  ## Generate auto paths to species shapefile
  my_sp <- as.character(Occs$species[1])
  
  pathSp <- paste("./data/rawdata/Biodiversity/iucn/", my_sp, sep = "", collapse = "")
  
  layerSp <- sub(".shp", "", list.files(path = pathSp, full.names = FALSE, pattern = ".shp")[1])
  
  
  ## Reading the shapefile of the species
  cat("Reading the IUCN shapefile\n")
  
  iucn <- readOGR(dsn = pathSp, layer = layerSp)
  
  ## Cropping to NC
  
  iucn <- crop(iucn, extent(eezNcPoly))
  
  ## Converting to raster
  cat("Converting shapefile to raster\n")
  
  rasIucnOccs <- raster(extent(iucn))
  
  res(rasIucnOccs) <- res(eezNcGrid)
  
  projection(rasIucnOccs) <- proj4string(eezNcGrid)
  
  origin(rasIucnOccs) <- origin(eezNcGrid)
  
  rasIucnOccs <- rasterize(iucn, field = "PRESENCE", rasIucnOccs) # assigning 1 for cells within the range
  
  rasIucnOccs <- reclassify(rasIucnOccs, cbind(1,2))
  
  
  cat("Merging eezNcGrid and IUCN rangemap raster\n")
  
  rasIucnOccs <- mask(rasIucnOccs, eezNcPoly)
  
  rasIucnOccs <- crop(rasIucnOccs, extent(eezNcPoly))
  
  antiAbs <- merge(rasIucnOccs,eezNcGrid)
  
  antiAbs <- reclassify(antiAbs, cbind(NA,2))
  
  antiAbs <- mask(antiAbs, eezNcPoly)
  
  antiAbs <- crop(antiAbs, extent(eezNcPoly))
  
  
  
  
  # Sampling the cells with the length of OccCells in the cells exterior to IUCN range and earthGrid
  
  OccsCells <- unique(Occs$cellNumber)
  
  cat("Sampling cells outside species range and earth\n")
  
  sampleClasses <- function(r, n)  {
    
    cellVal <- which(t(as.matrix(r)) == 1) # get All cells for class 2 (=Abs)
    
    samples <- sample(cellVal, n) # sample class's cell number
    
    return(samples)
  }
  
  
  pseudoAbsCells <- sampleClasses(antiAbs, length(OccsCells)) # multiply length(OccsCells) if we want another ratio than 50% of prevalence
  
  pseudoAbs <- as.data.frame(xyFromCell(antiAbs, pseudoAbsCells))
  
  colnames(pseudoAbs) <- c("decimalLongitude","decimalLatitude")
  
  pseudoAbs$cellNumber <- cellFromXY(eezNcGrid, as.matrix(pseudoAbs)) # Need to be corrected : lapply(pseudoAbs, cellFromXY(earthGrid, ) ... PseudoAbsCells arent the good ones
  
  pseudoAbs$species <- Occs$species[1]
  
  pseudoAbs$occurrence <- 0
  
  pseudoAbs <- pseudoAbs[,c(4,1,2,5,3)]
  
  pseudoAbs
  
  dataSet <- rbind(Occs, pseudoAbs) 
  
  cat("Done. Pseudo absences added to current dataset\n")
  
  dataSet
  
} # eo generateAbs


spDatasets <- lapply(spDatasets, generateAbs)


## binding substrate data


subsdata <- lapply(spDatasets, function(tabs, subst = rawSubs){
  
  subst <- subst[match(tabs$cellNumber, subst$cellNumber, nomatch=0),]
  
  tabs <- cbind(tabs, subst[ , -which(names(subst) %in% c("cellNumber"))])
  
}
)

testspt <- SpatialPoints(coords = as.matrix(spDatasets[[1]][,c(2,3)]), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") )




## Generating calib data for Travel_Dist


travelDist <- raster("./data/predictdata/Travel_Dist.tif")

travelDist <- as.data.frame(travelDist)
# travelDist <- cbind(coordinates(eezNcGrid), travelDist)
travelDist$cellNumber <- rownames(travelDist)


subsTravelData <- lapply(subsdata, function(tabs, trav = travelDist){
  
  trav <- trav[match(tabs$cellNumber, trav$cellNumber, nomatch=0),]
  
  trav$cellNumber = NULL
  
  tabs <- cbind(tabs, trav)
  
}
)


dir.create("./data/calibdata/regionmodel")

lapply(subsTravelData, function(tabs){
  
  write.csv(tabs, file = paste0("./data/calibdata/regionmodel/",tabs$species[1],"_speciesDataset_region.csv"), row.names = FALSE)
  
  return(paste0("dataset written in"," ","./data/regionmodel/bioclimodel/",tabs$species[1],"_speciesDataset_region.csv"))
  
})

test <- read.csv("./data/calibdata/regionmodel/Triaenodon_obesus_speciesDataset_region.csv")
# plot(travelDist)

