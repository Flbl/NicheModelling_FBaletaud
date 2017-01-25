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

OBISSP = function(my_sp, qc = c(1:7,10:19,21:30)){
  
my_sp = SpCheck(my_sp)

# Getting the occurrences

print("total number of records available")
my_sp_summ = checklist(my_sp)
print(my_sp_summ$records)

print("Getting data through Quality Flags :")
my_occs = occurrence(my_sp, qc = qc)

}


OBISSP("Carcharhinus amblyrhynchos")
OBISSP("Trianodon obesus")





SpCheck("Carcharhinus amblyrhynchos")
SpCheck("Carcharhinus amblyrhinchos")
SpCheck("Galeolamna tufiensis")
SpCheck("qsd")
my_sp <- "Carcharhinus amblyrhynchos"

SpCheck(my_sp)

my_sp_aphia <- get_wormsid(searchterm = my_sp, accepted = FALSE)
my_sp_taxo <- worms_records(ids = my_sp_aphia, marine_only = TRUE)
glimpse(my_sp_taxo) # Check status / valid_name / isMarine


# Getting occurences

my_occs_new = occurrence(my_sp, qc = c(1:7,10:19,21:30))
