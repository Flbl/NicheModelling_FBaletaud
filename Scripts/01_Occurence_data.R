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

# Merging the occurrences from OBIS and New Caledonia dataset

# Adding occurrence column to OBIS dataset
OBIS_occs$occurrence = 1

## Reading/formatting the dataset
NCRecords = read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE) 
NC_occs <- NCRecords[,c(2,3,7)]
NC_occs <- NC_occs[NC_occs$C_amblyrhynchos >= "1",]
colnames(NC_occs) <- c("decimalLatitude","decimalLongitude","individualCount")
NC_occs$occurrence = 1


# head(NC_occs)

## merge OBIS and NC

Occs <- merge(OBIS_occs, NC_occs, all = TRUE)

#Plot all the occurrences
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = Occs$decimalLatitude, lng = Occs$decimalLongitude), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)








# @knitr Manualchk

# Checking for same coordinates of points

Occs = distinct(Occs, decimalLongitude, decimalLatitude, .keep_all = TRUE)


# Double checking for points on earth using the bioclimatic grid

## Reading/creating the earth GRID
earthGrid<-nc_open("./environment/temp/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used


earthGrid<-ncvar_get(earthGrid,"thetao")
# temp<-as.data.frame(temp)

earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1
earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right

# plot(temp)

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus
writeRaster(earthGrid, "./environment/temp/earthGrid-1-12.asc",overwrite=T)

earthGrid = raster("./environment/temp/couche-1-12.asc")

## Gridding the occurrence dataset




## Multiplying the raster with the occurrence spatialpixelsdataframe




# Plotting for verification
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = Occs$decimalLatitude[1718], lng = Occs$decimalLongitude[1718]), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)


# @knitr datasetname

# Checking the datasets

table(Occs$datasetName)

# @knitr NTaiM

# Removing the National Taiwan Museum

Occs = Occs[-which(Occs$datasetName == "National Taiwan Museum"),]

nrow(Occs)



