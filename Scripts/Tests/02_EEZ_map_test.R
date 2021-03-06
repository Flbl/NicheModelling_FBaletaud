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
library(plyr)

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

# writeOGR(eezNcPoly, dsn = "./Miscellaneous", layer = "eezNcPoly", driver = "ESRI Shapefile")
# eezNcPoly <- readOGR(dsn = "./Miscellaneous", layer = "eezNcPoly")
# plot(earthGrid)
# plot(eezNcPoly, add = TRUE)





# intersect earthGrid with only the New Caledonian polygons

eezNcGrid <- mask(earthGrid, eezNcPoly)
eezNcGrid <- crop(eezNcGrid, extent(eezNcPoly))
eezNcGrid <- reclassify(eezNcGrid, cbind(1,NA)) # NA = earth/outlimits
eezNcGrid <- reclassify(eezNcGrid, cbind(2,1)) # 1 = Ocean
# writeRaster(eezNcGrid, "./Miscellaneous/eezNcGrid.asc", overwrite = TRUE)
# eezNcPoly <- readOGR(dsn = "./Miscellaneous", layer = "eezNcPoly")
# eezNcGridBehr <- eezNcGrid
# proj4string(eezNcGridBehr) <- "+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"


# eezNcSpGrid <- as(eezNcGrid, "SpatialGridDataFrame")
# proj4string(eezNcSpGrid) <- proj4string(eezNcGrid)

# eezNcPixGrid <- as(eezNcGrid, "SpatialPixelsDataFrame")
# proj4string(eezNcPixGrid) <- proj4string(eezNcGrid)


# plot(eezNcGrid)
# plot(eezNcPoly, add = TRUE)


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


# writeOGR(CoralNc, dsn = "./Miscellaneous", layer = "CoralNc", driver = "ESRI Shapefile")
# CoralNc <- readOGR(dsn = "./Miscellaneous", layer = "CoralNc")




# Merging the features of the polygons to one polygon

cropNcPoly <- function(shp, polyGrid = eezNcPolyGrid){
  
  shp <- gUnionCascaded(shp)
  
  shp <- gIntersection(eezNcPolyGrid, shp, byid = TRUE)
  
  shp
  
}

croppedNcPoly <- lapply(croppedNcshp, cropNcPoly)

croppedNcPoly <- croppedNcPoly[lapply(croppedNcPoly, length)!=0]

# Filtering the features through the occurrences and where they intersect (deleting features with no occurrences on them)
NCRecords <- read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE)

# keeping coords and species of interest
NCRecords <- NCRecords[,c(2:3, 7, 9, 11)]
NCRecords <- NCRecords[NCRecords$C_amblyrhynchos != 0 | NCRecords$C_melanopterus != 0 | NCRecords$T_obesus !=0,]
rownames(NCRecords) <- NULL

## creating spatial points for our occurrences

RecPoints <- SpatialPointsDataFrame(
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

RecPoly <- as(RecPoly, "SpatialPolygons")




# gIntersection(eezNcPolyGrid, RecPoints, byid = TRUE)


which(gIntersects(Canyons,RecPoints, byid = TRUE) == TRUE)


crossingFeatures <- lapply(croppedNcPoly, function(x, Rec = RecSpatial){ 
  
  true = which(gIntersects(x,RecSpatial, byid = TRUE) == TRUE)
  
  true
  
  })


crpdNcPolyFiltrd <- croppedNcPoly[lapply(crossingFeatures, length)!=0]



interCrss <- lapply(gIntersection(crpdNcPolyFiltrd,eezNcPolyGrid, byid = TRUE))


# 167.93750001999999 -21.354166670000001 at 167.93750001999999 -21.354166670000001 
errdat <- data.frame(lat = c(167.93750001999999,167.93750001999999 ), long =c(-21.354166670000001, -21.354166670000001))

errCoords <- coordinates(errdat)
errSpat <- SpatialPoints(coords = errCoords,
                         proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



InterCanyons


Cells <- extract(eezNcGrid, InterCanyons, cellnumbers = TRUE)



# Extraire le labpt (centroid) de chaque features de l'objet InterCanyons et faire un extract dessus pour renvoyer la cellule


coordinates(InterCanyons[1:5])

InterCanyons@polygons[[5]]@Polygons[[1]]@labpt

dim(getSpPPolygonsLabptSlots(InterCanyons))

Cells <- extract(eezNcGrid, coordinates(InterCanyons), cellnumber = TRUE)
dim(Cells)






# extract percentages of each polygons in each cells





######################### Bordel ######################################
test <- area(eezNcGrid, na.rm = TRUE)


# Change the projection of everything to calculate areas

test <- interCrss[[1]]

test <- spTransform(test, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))

# test <- spTransform(eezNcPoly, CRS("+init=epsg:2984"))

gArea(test, byid = TRUE)/(1000*1000)

plot(NcGridCea)


# Creating the final data.frame with cells numbers and values associated

gridCells <- as.data.frame(eezNcGrid)

abyssC <- extract(eezNcGrid, coordinates(interCrss[[1]]), cellnumber = TRUE)
abyssC <- as.data.frame(abyssC)


# extract percentages of each polygons in each cells

test <- spTransform(test, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))

# test <- spTransform(eezNcPoly, CRS("+init=epsg:2984"))

berhSubst <- lapply(interCrss, spTransform, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")) 

abyssA <- gArea(berhSubst[[1]], byid = TRUE)/(1000*1000)

names(abyssA) <- abyssC[,"cells"]


abyssC <- extract(eezNcGrid, coordinates(interCrss[[1]]), cellnumber = TRUE)
abyssC <- as.data.frame(abyssC)

behrSubst <- lapply(interCrss, spTransform, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs")) 

abyss[!is.na(featCov[[1]][,1])] <- featCov[[1]][,3]

abyss <- replace(abyss, abyss[which(names(abyss) == featCov[[1]][,1])] ,  featCov[[1]][,3])


featCovRep <- featCov[[1]][,1]
length(featCovRep) <- length(abyss)
featCovRep[is.na(featCovRep)] <- 0

abyss[abyss == featCov[[1]][,1]] <- featCov[[1]][,3]


abyss <- replace(abyss, abyss[which(names(abyss) == featCov[[1]][,1])] ,  featCov[[1]][,3])



names(abyss) <- 1:dim(substrateNoNa)[1]

abyss <- replace(abyss, match(featCov[[1]][,1], names(abyss)),  featCov[[1]][,3])

mtch <- match(featCov[[1]][,1], names(abyss))

abyss[featCov[[1]][,1] %in% names(abyss)] <- featCov[[1]][,3]


featCov[[1]][,1] %in% names(abyss)

match(featCov[[1]][,1], names(abyss))

# abyss <- featCov[[1]][,3]

abyss <- data.frame(area = rep (0, dim(substrateNoNa)[1]))

abyss$area[which(rownames(abyss) == match(featCov[[1]][,1], rownames(abyss)))] <- featCov[[1]][,3]

mtch <- as.vector(match(featCov[[1]][,1], rownames(abyss)))

abyss[which(abyss$cells == mtch)] <- featCov[[1]][,3]

abyss$cells <- rownames(abyss)

abyss = merge(abyss, featCov[[1]][,c(1,3)], by = )# match(featCov[[1]][,1], rownames(abyss)))

abyss = cbind(substrateNoNa, abyss)




abyss = cbind.fill(abyss, featCov[[1]][,c(1,3)])
length(abyss) <- nrow(substrateNoNa)
abyss <- merge(abyss, featCov[[1]][,c(1,3)], by = intersect(rownames(abyss),featCov[[1]][,1] ) )

# names(abyss) = c(1:length(abyss))

# abyss$rownames <- rownames(abyss)

abyss[abyss$rownames==featCells[[1]][,"cells"],]

abyss <- replace(abyss, names(abyss)==featCells[[1]][,"cells"], featCov[[1]][,3] )


rownames(abyss) <- featCells[[1]][,"cells"]

abyss[rownames(featCov[[1]])] <- featCov[[1]][,3]


abyssCov <- cbind(abyssC, abyssA) 

abyssCov1 <- rep (NA, dim(gridCells)[1])

abyssCov1[!is.na(gridCells)] <- abyssCov[,3]

abyssCov <- cbind(gridCells, abyssCov1)




featCov1 <- lapply(featCov, function(x){
  
  y <- names(featCov)
  
  # res <- grep(,)


x <- rep (NA, dim(gridCells)[1])

x[!is.na(gridCells)] <- featCov[[y]][,"area"] == TRUE

x

})


featCov1 <- for(i in length(featCov)){
  
  
  featCov1[[i]][!is.na(gridCells)] <- featCov[[i]][,"area"]
  
}










ncGridValues <- extract(eezNcGrid, croppedNcPoly[[4]], small = TRUE)


## Get class counts for each polygons
valCounts <- lapply(ncGridValues, table)

## Calculate class percentage

valPct <- lapply(valCounts, FUN=function(x, eezNcGrid){ x / sum() } )


#########################################################


getSingleGridData <- function(x, gridd = gridCells){
  
  substrate <- gridd
  
  substrateNoNa <- data.frame(layer = gridd$layer[!is.na(gridd$layer)])
  
  rownames(substrateNoNa) <- na.omit(rownames(as.data.frame(substrate))[substrate[["layer"]] == 1])
  
  substrateNoNa$RN <- as.numeric(rownames(substrateNoNa))
  
  abyss1 <- data.frame(area = x[,3], RN = x[,1])
  
  abyss <- aggregate(abyss1, by = list(cells = abyss1$RN), sum)
  
  rownames(abyss) <- abyss$cells
  
  zeroNames <- substrateNoNa$RN[is.na(match(substrateNoNa$RN, abyss$cells))]
  
  zeroDat <- data.frame(area = numeric(dim(substrateNoNa)[1] - dim(abyss)[1]), row.names = zeroNames)
  
  abyss <- data.frame(area = abyss$area, row.names = rownames(abyss))
  
  abyss <- rbind(abyss, zeroDat)
  
  abyss <- data.frame(area = abyss[ order(as.numeric(row.names(abyss))), ])
  
  rownames(abyss) <- rownames(substrateNoNa)
  
  abyssNa <- rep (NA, dim(substrateNoNa)[1])
  
  abyssNa[!is.na(gridd)] <- abyss[,"area"]
  
  abyssNa
  
}#eo getSingleGridData

# substData <- cbind(substrate, abyssNa)

# substData
x = featCov[[2]]
test <- getSingleGridData(featCov[[2]])

test <- lapply(featCov,getSingleGridData)







testSubs <- cbind(coordinates(eezNcGrid), substrate)
testSubs <- na.omit(testSubs)

testSpt <-  SpatialPointsDataFrame(
  coords = cbind(testSubs$x, testSubs$y),
  data = testSubs,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


testRast <- raster(extent(eezNcGrid))
res(testRast) <- res(eezNcGrid)
projection(testRast) <- proj4string(eezNcGrid)
origin(testRast) <- origin(eezNcGrid)

testRast <- rasterize(testSpt,field = "Ridges", testRast)
plot(testRast)

testRast <- raster(test)
test <- as.matrix(substrate$Abyss)





substrateSpatial <-  SpatialPointsDataFrame(
  coords = cbind(substrate$x, substrate$y),
  data = substrate,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


featRast <- raster(extent(eezNcGrid))
res(featRast) <- res(eezNcGrid)
projection(featRast) <- proj4string(eezNcGrid)
origin(featRast) <- origin(eezNcGrid)

RidgesRast <- rasterize(substrateSpatial ,field = "Ridges", featRast)
GuyotsRast <- rasterize(substrateSpatial ,field = "Guyots", featRast)
CoralRast <- rasterize(substrateSpatial ,field = "Coral", featRast)
SlopeRast <- rasterize(substrateSpatial ,field = "Slope", featRast)
AbyssRast <- rasterize(substrateSpatial ,field = "Abyss", featRast)
ShelfRast <- rasterize(substrateSpatial ,field = "Shelf", featRast)
ShelfRast <- rasterize(substrateSpatial ,field = "Escarpments", featRast)



par(mfrow = c(2,2))
plot(RidgesRast)
plot(GuyotsRast)
plot(CoralRast)
plot(SlopeRast)
plot(AbyssRast)
plot(ShelfRast)






plot(eezNcPolyGrid)

lapply(croppedNcPoly, plot, col = sample(colors(), 1), add = TRUE)

# plot(croppedNcPoly[[1]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[2]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[3]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[4]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[5]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[6]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[7]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[8]], add = TRUE, col = sample(colors(), 1))
# plot(croppedNcPoly[[10]], add = TRUE, col = sample(colors(), 1))


# par(mfrow = c(3,2))
# par(oma = c(0.1, 0.1, 0.1, 0.1))
layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = TRUE), widths = 5, heights = 2 )












# test 


InterCanyons = readOGR(dsn = "./Environment/geomorph", layer = "Canyons")

InterCanyons = crop(InterCanyons, extent(eezNcPolyGrid))

InterCanyons

Canyons = gUnionCascaded(InterCanyons)

InterCanyons <- gIntersection(Canyons,eezNcPolyGrid, byid = TRUE)

InterCanyons


Cells <- extract(eezNcGrid, InterCanyons, cellnumbers = TRUE)

InterCanyons@polygons[[1]]@Polygons[[1]]@coords

# test <- SpatialPolygons(InterCanyons@polygons[[1]]@Polygons)



# test 2


poly = crpdNcPolyFiltrd[[3]]

poly <- gIntersection(poly,eezNcPolyGrid, byid = TRUE)

plot(eezNcGrid)
plot(poly, add = TRUE, col = "red")





polytest1coords <- InterCanyons@polygons[[1]]@Polygons[[1]]@coords

polyTest1 <- SpatialPointsDataFrame(
  coords = coordinates(polytest1coords),
  data = data.frame(polytest1coords),
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


writeOGR(polyTest1, dsn = "./Miscellaneous", layer = "polyTest1", driver = "ESRI Shapefile")
writeOGR(eezNcPolyGrid, dsn = "./Miscellaneous", layer = "eezNcPolyGrid", driver = "ESRI Shapefile")



length(extract(eezNcGrid, InterCanyons, small = TRUE))


Cells[[996]]
InterCanyons@polygons[[996]]@Polygons[[1]]@coords

eezNcGridTest <- eezNcGrid

eezNcGridTest[20249] <- 2
eezNcGridTest[20250] <- 3
eezNcGridTest[20466] <- 4
eezNcGridTest[20465] <- 5

plot(eezNcGridTest)
plot(InterCanyons, add = TRUE)
polytest996coords <- InterCanyons@polygons[[996]]@Polygons[[1]]@coords

polyTest996 <- SpatialPointsDataFrame(
  coords = coordinates(polytest996coords),
  data = data.frame(polytest996coords),
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


writeOGR(polyTest996, dsn = "./Miscellaneous", layer = "polyTest996", driver = "ESRI Shapefile")


Cells[[351]]


polytest351coords <- InterCanyons@polygons[[351]]@Polygons[[1]]@coords


polyTest351 <- SpatialPointsDataFrame(
  coords = coordinates(polytest351coords),
  data = data.frame(polytest351coords),
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

writeOGR(polyTest351, dsn = "./Miscellaneous", layer = "polyTest351", driver = "ESRI Shapefile")

eezNcGridTest[13229] <- 2
eezNcGridTest[13228] <- 3
eezNcGridTest[13227] <- 4

plot(eezNcGridTest)
plot(InterCanyons, add = TRUE)

centro <- gCentroid(InterCanyons)
plot(centro, col = "red", add = TRUE, cex = 5)


InterCanyons@polygons[[351]]@Polygons[[1]]@labpt

# gArea(InterCanyons)
# 
# Canyons
# plot(Canyons)
# 
# 
# vals <- extract(eezNcGrid, Canyons)




r <- raster(ncol=36, nrow=18)
r[] <- round(runif(ncell(r),1,10),digits=0)


cds1 <- rbind(c(-180,-20), c(-160,5), c(-60, 0), c(-160,-60), c(-180,-20))
cds2 <- rbind(c(80,0), c(100,60), c(120,0), c(120,-55), c(80,0))
polys <- SpatialPolygonsDataFrame(SpatialPolygons(list(Polygons(list(Polygon(cds1)), 1),
                                                       Polygons(list(Polygon(cds2)), 2))),data.frame(ID=c(1,2)))

v <- extract(r, polys, small = TRUE)

v.counts <- lapply(v,table)

v.pct <- lapply(v.counts, FUN=function(x){ x / sum(x) } )

class.df <- as.data.frame(t(sapply(v.pct,'[',1:length(unique(r)))))



# Intersecting the variables polygons to the eezNcPolyGrid

gArea(test, croppedNcshp[[4]], byid = TRUE)

gIntersection(eezNcPolyGrid, croppedNcshp[[4]], byid = TRUE)

InterCanyons <- gIntersection(eezNcPolyGrid, croppedNcshp[[4]], byid = TRUE)

plot(InterCanyons, add = TRUE, col = "red")











# eezCoralClip <- gIntersection(eezNcGrid, CoralNc, byid = TRUE) # Ne fonctionne pas en raster,polygons





# Abyss <- cropNc("Abyss")
# 
# Abyssal_Classification <- cropNc("Abyssal_Classification")
# 

# Shelf_Classification <- cropNc("Shelf_Classification")

# 
# cropNc("Seamounts")
# 
# plot(Abyss)

##################################################

# Whole 04_EEZ_map.R script before modifying it to get all selected features instead
# of filtering through features intersecting with occurrences

##################################################






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

croppedNcshp <- croppedNcshp[lapply(croppedNcshp, length)!=0]


#Adding the coral

## Crop the Coral 

CoralShp = readOGR("./data/rawdata/Environment/WCMC008_CoralReef2010_v1_3/01_Data", layer = "14_001_WCMC008_CoralReef2010_v1_3")

CoralNc <- crop(CoralShp, extent(eezNcPoly))


croppedNcshp$Coral <- CoralNc


# Merging the features of the polygons to one polygon # Changed due to issues in the resulting polygons

croppedNcPoly <- croppedNcshp # before lapply(croppedNcshp, gUnionCascaded)

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





