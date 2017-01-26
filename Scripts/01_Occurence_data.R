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

my_occs = occurrence(my_sp, qc = qc)

assign("my_occs", my_occs, .GlobalEnv)

}



OBISSP("Carcharhinus amblyrhynchos")
# OBISSP("Trianodon obesus")

# Plotting
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = my_occs$decimalLatitude, lng = my_occs$decimalLongitude), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)


# @knitr Manualchk


# Removing points not in the oceans (using ecoregions shapefile)

Ecoregion = readOGR(dsn = "./Environment/ecoregion/terEcorDissolved", layer = "terEcorDissolved")

my_occs_spatial <- SpatialPointsDataFrame(
  coords = cbind(my_occs$decimalLongitude, my_occs$decimalLatitude),
  data = my_occs,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

is_occ_Marine <- !(gIntersects(my_occs_spatial,Ecoregion, byid = TRUE))

my_occs = my_occs[is_occ_Marine == TRUE, ]

# Plotting for verification
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = my_occs$decimalLatitude, lng = my_occs$decimalLongitude), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)


# Checking the datasets

table(my_occs$datasetName)

# Removing the National Taiwan Museum

my_occs = my_occs[-which(my_occs$datasetName == "National Taiwan Museum"),]

nrow(my_occs)

# @knitr mergeOccs

# Merging the occurrences from OBIS and New Caledonia dataset

# Adding occurrence column to OBIS dataset
my_occs$occurrence = 1

## Reading/formatting the dataset
NCRecords = read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE) 
NC_occs <- NCRecords[,c(2,3,7)]
NC_occs <- NC_occs[NC_occs$C_amblyrhynchos == "1",]
colnames(NC_occs) <- c("decimalLatitude","decimalLongitude","occurrence")

head(NC_occs)

## merge OBIS and NC

my_occs_merged <- merge(my_occs, NC_occs, all = TRUE)

#Plot all the occurrences
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = my_occs_merged$decimalLatitude, lng = my_occs_merged$decimalLongitude), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)

