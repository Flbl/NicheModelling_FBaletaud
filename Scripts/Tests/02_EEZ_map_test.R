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


err <- 
InterCanyons


Cells <- extract(eezNcGrid, InterCanyons, cellnumbers = TRUE)






# extract percentages of each polygons in each cells

ncGridValues <- extract(eezNcGrid, croppedNcPoly[[4]], small = TRUE)


## Get class counts for each polygons
valCounts <- lapply(ncGridValues, table)

## Calculate class percentage

valPct <- lapply(valCounts, FUN=function(x, eezNcGrid){ x / sum() } )


















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


# Extraire le labpt (centroid) de chaque features de l'objet InterCanyons et faire un extract dessus pour renvoyer la cellule


coordinates(InterCanyons[1:5])

InterCanyons@polygons[[5]]@Polygons[[1]]@labpt

dim(getSpPPolygonsLabptSlots(InterCanyons))

Cells <- extract(eezNcGrid, coordinates(InterCanyons), cellnumber = TRUE)
dim(Cells)



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

gArea(eezNcPolyGrid, croppedNcshp[[4]], byid = TRUE)

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







