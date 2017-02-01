#################################################################
##                                                             ##
##            Getting and formatting occurence data            ##         
##                                                             ##
#################################################################

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


# knitr NCdata



# @knitr OBIS

############  OBIS  #############


# Checking the name
SpCheck<- function(checkSp){
  checkSp_aphia <- get_wormsid(searchterm = checkSp, accepted = FALSE)
  checkSp_taxo <- worms_records(ids = checkSp_aphia, marine_only = TRUE)
  
  if(is.na(checkSp_aphia)){
    stop("Species not found")
  }
  
  if(checkSp_taxo$status == "unaccepted"){
    print("invalid name : assigning latest name")
    my_sp = checkSp_taxo$valid_name
  }
  if(checkSp_taxo$status == "accepted"){
    print("Valid name")
    my_sp = checkSp
    
  }
  review = data.frame(Entered_Name = checkSp, Reason = checkSp_taxo$unacceptreason, Accepted_Name = my_sp)
  print(review)
  return(my_sp)
  
}

# SpCheck("Carcharhinus amblyrhynchos")
# SpCheck("Carcharhinus amblyrhinchos")
# SpCheck("Galeolamna tufiensis")
# SpCheck("Internecivus raptus") 

OBISSP = function(my_sp, qc = c(1:7,10:19,21:30)){
  
my_sp = SpCheck(my_sp)

# Getting the occurrences

print("total number of records available")
my_sp_summ = checklist(my_sp)
print(my_sp_summ$records)

print("Getting data through Quality Flags :")

OBIS_occs = occurrence(my_sp, qc = qc)

assign("OBIS_occs", OBIS_occs, .GlobalEnv)

}


OBISSP("Carcharhinus amblyrhynchos")
# OBISSP("Trianodon obesus")

#######################################################


# @knitr plotocc1

# Plotting
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = OBIS_occs$decimalLatitude, lng = OBIS_occs$decimalLongitude), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)


#######################################################


# @knitr mergeOccs

mergeOccs = function(my_sp){

# Calling the OBISSP function
  
  OBISSP(my_sp)
  
# Merging the occurrences from OBIS and New Caledonia dataset

# Adding occurrence column to OBIS dataset
OBIS_occs$occurrence = 1

## Reading/formatting the dataset
NCRecords <- read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE) 
match = paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2], sep = "", collapse = "")
matchCol = grep(match,names(NCRecords))
names(NCRecords)[matchCol] <- my_sp

NC_occs <- NCRecords[,c(2,3, matchCol)]
NC_occs <- NC_occs[NC_occs[,3] >= "1",]
colnames(NC_occs) <- c("decimalLatitude","decimalLongitude","individualCount")
NC_occs$occurrence = 1
NC_occs$datasetName <- "APEX"
rownames(NC_occs) <- NULL # Reassign the rownames of NC_occs


assign("NC_occs", NC_occs, .GlobalEnv)


## merge OBIS and NC

OccsRaw <- merge(OBIS_occs, NC_occs, all = TRUE)
assign("OccsRaw", OccsRaw, .GlobalEnv)

Occs <- OccsRaw[, names(NC_occs)]
assign("Occs", Occs, .GlobalEnv)


}

mergeOccs("Carcharhinus amblyrhynchos")

#Plot all the occurrences
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = Occs$decimalLatitude, lng = Occs$decimalLongitude), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)








# @knitr Manualchk

Manualchk = function(my_sp){
  
  mergeOccs(my_sp)


# Checking for same coordinates of points

Occs = distinct(Occs, decimalLongitude, decimalLatitude, .keep_all = TRUE)


# checking for points on earth using the bioclimatic grid

## Reading/creating the earth GRID
earthGrid<-nc_open("./environment/temp/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used


earthGrid<-ncvar_get(earthGrid,"thetao")
# temp<-as.data.frame(temp)


earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1
earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus
# writeRaster(earthGrid, "./environment/temp/earthGrid-1-12.asc",overwrite=T)
# earthGrid = raster("./environment/temp/earthGrid-1-12.asc")
# plot(earthGrid)


## creating spatial points for our occurrences

occsSpatial = SpatialPointsDataFrame(
  coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
  data = Occs,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)

# plot(occsSpatial, add = TRUE)

Occs <- Occs[-which(raster::extract(earthGrid,occsSpatial) == 1),]

# writeOGR(occsSpatial, dsn = "./environment/temp", layer = "OccsSpatial302", driver = "ESRI Shapefile")


rownames(Occs) <- NULL # Reassign the rownames of Occs


# Double checking with the GSHHS "special islets" regions shapefile

occsSpatial = SpatialPointsDataFrame(
  coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
  data = Occs,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)

earthRegion <- readOGR(dsn = "./Environment/GSHHS_region", layer = "GSHHS_f_L1")

is_occ_terrest <- gIntersects(occsSpatial,earthRegion, byid = TRUE, returnDense = FALSE)



occ_terrest <- as.numeric(names(unlist(is_occ_terrest)))

Occs <- Occs[-occ_terrest,]

assign(my_sp, Occs, .GlobalEnv)

}


Manualchk("Carcharhinus amblyrhynchos")









GetCleanedOcc = function(my_sp, qc = c(1:7,10:19,21:30)){
  
  # CheckSp function
  checkSp <- gsub("_", " ",my_sp)
  
  
  checkSp_aphia <- get_wormsid(searchterm = checkSp, accepted = FALSE)
  checkSp_taxo <- worms_records(ids = checkSp_aphia, marine_only = TRUE)
  
  if(is.na(checkSp_aphia)){
    stop("Species not found")
  }
  
  if(checkSp_taxo$status == "unaccepted"){
    print("invalid name : assigning latest name")
    my_sp <- checkSp_taxo$valid_name
  }
  if(checkSp_taxo$status == "accepted"){
    print("Valid name")
    my_sp <- checkSp
    
  }
  review = data.frame(Entered_Name = checkSp, Reason = checkSp_taxo$unacceptreason, Accepted_Name = my_sp)
  print(review)
  
  
  
  # OBISSP function
  
  print("total number of records available")
  my_sp_summ = checklist(my_sp)
  print(my_sp_summ$records)
  
  print("Getting data through Quality Flags :")
  
  OBIS_occs = occurrence(my_sp, qc = qc)
  
  print("Keeping only occurrences after 1980")
  OBIS_occs = OBIS_occs[OBIS_occs$yearcollected >= 1980,]
  
  if (nrow(OBIS_occs) <= 2) stop("Not enough records")
  
  #assign(paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2],"OBIS_occs", sep = "", collapse = ""), OBIS_occs, .GlobalEnv)
  
  
  
  # Merging the occurrences from OBIS and New Caledonian dataset
  
  ## Adding occurrence column to OBIS dataset
  OBIS_occs$occurrence = 1
  
  ## Reading/formatting the dataset
  NCRecords <- read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE) 
  print("Creating/Formatting the NC dataset")
  match = paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2], sep = "", collapse = "")
  matchCol = grep(match,names(NCRecords))
  names(NCRecords)[matchCol] <- my_sp
  
  NC_occs <- NCRecords[,c(2,3, matchCol)]
  NC_occs <- NC_occs[NC_occs[,3] >= "1",]
  colnames(NC_occs) <- c("decimalLatitude","decimalLongitude","individualCount")
  NC_occs$occurrence = 1
  NC_occs$datasetName <- "APEX"
  rownames(NC_occs) <- NULL # Reassign the rownames of NC_occs
  
  
  #assign(paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2],"NC_occs", sep = "", collapse = ""), NC_occs, .GlobalEnv)
  
  
  ## merge OBIS and NC
  
  print("Merging datasets")
  OccsRaw <- merge(OBIS_occs, NC_occs, all = TRUE)

  Occs <- OccsRaw[, names(NC_occs)]

  
  # Manual check of the dataset
  
  print("deleting coordinates duplicates")
  
  if(any(is.na(Occs$decimalLatitude) & is.na(Occs$decimalLongitude))) Occs = Occs[-which(is.na(Occs$decimalLatitude) & is.na(Occs$decimalLongitude)),]
  
  Occs = distinct(Occs, decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  
  
  print("replacing NA with 1 to rows without individual count (fisheries)")
  
  Occs[is.na(Occs$individualCount),3] <- 1
  
  
  
  rownames(Occs) <- NULL # Reassign the rownames of NC_occs
  
  # checking for points on earth using the bioclimatic grid
  
  print("Check points on earth using bioclimatic grid")
  
  # ## Reading/creating the earth GRID
  # earthGrid<-nc_open("./environment/temp/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used
  # 
  # 
  # earthGrid<-ncvar_get(earthGrid,"thetao")
  # # temp<-as.data.frame(temp)
  # 
  # 
  # earthGrid[(is.na(earthGrid))]<-1000
  # earthGrid[(earthGrid<1000)]<-NA
  # earthGrid[(earthGrid == 1000)]<-1
  # earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right
  # 
  # earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus
  # writeRaster(earthGrid, "./environment/temp/earthGrid-1-12.asc",overwrite=T)
  # earthGrid = raster("./environment/temp/earthGrid-1-12.asc")
  # plot(earthGrid)
  
  
  ## creating spatial points for our occurrences
  
  occsSpatial = SpatialPointsDataFrame(
    coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
    data = Occs,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  )
  
  # plot(occsSpatial, add = TRUE)
  
  Occs <- Occs[-which(raster::extract(earthGrid,occsSpatial) == 1),]
  
  # writeOGR(occsSpatial, dsn = "./environment/temp", layer = "OccsSpatialGcuvier", driver = "ESRI Shapefile")
  
  
  rownames(Occs) <- NULL # Reassign the rownames of Occs
  
  
  
  
  
  # Double checking with the GSHHS "special islets" regions shapefile
  
  print("Check points on earth using GSHHS shapefile")
  
  
  occsSpatial = SpatialPointsDataFrame(
    coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
    data = Occs,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  )
  
  # earthRegion <- readOGR(dsn = "./Environment/GSHHS_region", layer = "GSHHS_f_L1")
  
  is_occ_terrest <- gIntersects(occsSpatial,earthRegion, byid = TRUE, returnDense = FALSE)
  
  occ_terrest <- as.numeric(names(unlist(is_occ_terrest)))
  
  Occs <- Occs[-occ_terrest,]
  
  # Adding the species name
  Occs$species <- gsub(" ", "_",my_sp)
  
  #my_sp_under <- gsub(" ", "_",my_sp)
  #assign(my_sp_under, Occs, .GlobalEnv)
  
  Occs
  
}#eo GetCleanedOcc







## Reading/creating the earth GRID

earthGrid<-nc_open("./environment/temp/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used


earthGrid<-ncvar_get(earthGrid,"thetao")



earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1
earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus

earthRegion <- readOGR(dsn = "./Environment/GSHHS_region", layer = "GSHHS_f_L1")





test <- GetCleanedOcc("Galeocerdo_cuvier")

ensemblecheck <- function(species){
  
  check <- lapply(species, checklist, qc = c(1:7,10:19,21:30))
  
  names(check) <- species

  check <- do.call(rbind,check) # -> rbind tous les tableaux de la liste
  
  check = data.frame(species = check$species, records = check$records)
  
  
}

species <- c("Carcharhinus amblyrhynchos","Galeocerdo cuvier","Nebrius ferrugineus","Sphyrna lewini")


cleandOccs <- lapply(species,GetCleanedOcc)

names(cleandOccs) <- species

cleandOccs[["Carcharhinus_amblyrhynchos"]] # <=> cleandOccs$Carcharhinus_amblyrhynchos

cleandOccsMerged <- do.call(rbind,cleandOccs) # -> rbind tous les tableaux de la liste

# GetCleanedOcc("Carcharhinus_amblyrhynchos")
# 
# GetCleanedOcc("Galeocerdo cuvier")
# 
# GetCleanedOcc("Nebrius ferrugineus")
# 
# GetCleanedOcc("Sphyrna lewini")









plot(earthGrid)
plot(occsSpatial, add = TRUE)
plot(occsSpatial[occ_terrest,], add = TRUE, col = "red")


# Plotting for verification
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = Occs$decimalLatitude[21], lng = Occs$decimalLongitude[21]), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)

# is_occ_terrest_ids <- lapply(is_occ_terrest,function(x){
#   
#   cat(x,"\n")
#   
#   if(is.null(x)) return(NULL)
#   
#   id <- ((earthRegion@polygons)[[x]])@ID
#   
#   return(id)
#   
#   
# })



my_sp_under <- gsub(" ", "_",my_sp)

