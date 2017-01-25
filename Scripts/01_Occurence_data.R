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

# @knitr GBIF

############  GBIF  #############

ambly = name_lookup("Carcharhinus amblyrhynchos", rank = "species", return = "data")


name_suggest("Carcharhinus amblyrhynchos", rank = "species") #$key[1] for first suggested key (not necessarly the good one)
occ_search(scientificName = "Carcharhinus amblyrhynchos", limit = 20)
occ_count(taxonKey = 2418064, georeferenced = TRUE)

amblyGbif = occ_search(scientificName = "Carcharhinus amblyrhynchos", return = "all")
amblyGbif = occ_search(scientificName = "Carcharhinus amblyrhynchos", return = "meta")

amblyGbif = occ_search(scientificName = "Carcharhinus amblyrhynchos", return = "data", field = c("name","decimalLongitude","decimalLatitude","basisOfRecord","coordinateUncertaintyInMeters","depth","depthAccuracy","waterBody","year","locality","samplingProtocol"))

gbifmap(amblyGbif)

#################################


# @knitr OBIS

############  OBIS  #############

# Checking the name

my_sp <- "Carcharhinus amblyrhynchos"
my_sp_aphia <- get_wormsid(searchterm = my_sp, accepted = FALSE)
my_sp_taxo <- worms_records(ids = my_sp_aphia, marine_only = TRUE)
glimpse(my_sp_taxo) # Check status / valid_name / isMarine
my_sp_summ = checklist(my_sp)
my_sp_summ$records

# Getting occurences

my_occs <- occurrence(my_sp)
bb_occs <- bbox(cbind(my_occs$decimalLongitude, my_occs$decimalLatitude))
bb_occs

world = map_data("world")

worldmap <- ggplot(world, aes(x=long, y=lat)) +
  geom_polygon(aes(group=group)) +
  scale_y_continuous(breaks = (-2:2) * 30) +
  scale_x_continuous(breaks = (-4:4) * 45) +
  theme(panel.background = element_rect(fill = "steelblue")) +
  coord_equal()

# Get the centroid of the whole coordinates
centro <- data.frame(Long = mean(my_occs$decimalLongitude), Lat = mean(my_occs$decimalLatitude))

occ_map = worldmap + geom_point(data = my_occs, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkorange", shape = 21, alpha = 2/3) + geom_point(data = centro, aes(x = Long, y = Lat), colour = "red", shape = 21, alpha = 2/3)

occ_map
# occ_map + coord_map("sinusoidal")
table(my_occs$geodeticDatum)

# Basemap from google maps (better for regional maps)
my_map <- get_map(location = bb_occs, maptype = "satellite")
ggmap(my_map) + geom_point(data = my_occs, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkorange", shape = 21, alpha = 2/3)

# Plot with a colour code by the original scientific name recorded (OBIS translated it into the accepted name)

worldmap + geom_point(data = my_occs, aes(
  x = decimalLongitude, 
  y = decimalLatitude, 
  colour = originalScientificName), 
  shape = 21, 
  alpha = 2/3)

# Plot with decade collected
my_occs$decade <- with(my_occs, 10*round(yearcollected/10, 0))
worldmap + geom_point(data = my_occs,
                      aes(x = decimalLongitude, y = decimalLatitude, colour = decade), 
                      shape = 21, 
                      alpha = 2/3) + scale_colour_gradient(low = "white", high = "darkorange")

# In case of numerous occurences for a species :
# occurence(my_sp, field = c("species", "decimalLongitude", "decimalLatitude"))

# Quality control

glimpse(my_occs)
# Getting an idea of the proportion of records with data for each field

round(sort(apply(!is.na(my_occs), 2, mean)), 3)

table(my_occs$originalScientificName)

# Filter by QCflags

# See which flags are on
as.logical(intToBits(my_occs$qc[1]))

filter_by_QcFlags <- function(occ_dat, qc_var = "qc", qc_flags){
  get_allon_ids <- function(qc_var, qc_flags) {
    
    mask <- 2^(qc_flags - 1)
    qc_flags_on <- sapply(qc_var,function(x){sum(bitwAnd(x,mask) > 0)})
    all_on <- which(qc_flags_on == length(qc_flags))
    all_on
  } 
  
  if(min(qc_flags, na.rm = TRUE) < 1 | max(qc_flags, na.rm = TRUE) >30 | !(class(qc_flags) %in% c("numeric","integer"))){
    stop("Invalid values for qc_flags, must be integers in the range 1:30", call. = FALSE)
  }
  
  
  if(min(qc_flags, na.rm = T) < 1 | max(qc_flags, na.rm = T) > 30 |
     !(class(qc_flags) %in% c("numeric", "integer"))){
    stop("Invalid values for qc_flags, must be integers in the range 1:30",
         call. = FALSE)
  }
  
  
  if(sum(c(8, 9, 20) %in% qc_flags) > 0){
    stop("Flags 7,8,9,16,20,25,26,27,28 are currently disabled and no records would be returned by your query",
         call. = FALSE)
  }
  
  if(qc_var != "qc"){occ_dat <- plyr::rename(occ_dat, setNames('qc', eval(qc_var)))}
  id_all_on <- get_allon_ids(occ_dat$qc, qc_flags)
  
  occ_dat <- occ_dat[id_all_on, ]
  
  if(qc_var != "qc"){occ_dat <- plyr::rename(occ_dat, setNames(eval(qc_var), 'qc'))}
  
  return(occ_dat)
}

my_occs_filt = filter_by_QcFlags(my_occs, qc_flags = c(1:7,10:19,21:30))
nrow(my_occs) - nrow(my_occs_filt)

centrofilt = data.frame(Long = mean(my_occs_filt$decimalLongitude), Lat = mean(my_occs_filt$decimalLatitude))

  
occ_map_filt = worldmap + geom_point(data = my_occs_filt, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkorange", shape = 21, alpha = 2/3) + geom_point(data = centrofilt, aes(x = Long, y = Lat), colour = "red", shape = 21, alpha = 2/3)
occ_map_filt

# my_occs$qcnum = qcflags(my_occs$qc,c(1:7,10:19,21:30))
# colors <- c("#ee3300", "#86b300")[my_occs$qcnum + 1]
# 
# leaflet() %>%
#   addProviderTiles("CartoDB.Positron") %>%
#   addCircleMarkers(popup = paste0(my_occs$datasetName, "<br/>", my_occs$catalogNumber, "<br/><a href=\"http://beta.iobis.org/explore/#/dataset/", my_occs$resourceID, "\">OBIS dataset page</a>"), data = data.frame(lat = my_occs$decimalLatitude, lng = my_occs$decimalLongitude), radius = 3.5, weight = 0, fillColor = colors, fillOpacity = 1)

