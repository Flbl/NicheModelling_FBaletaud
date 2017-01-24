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

occ_map = worldmap + geom_point(data = my_occs, aes(x = decimalLongitude, y = decimalLatitude), colour = "darkorange", shape = 21, alpha = 2/3)

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


amblyObis$qcnum <- qcflags(amblyObis$qc, c(28))
colors <- c("#ee3300", "#86b300")[sturgeon_data$qcnum + 1]

leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(popup = paste0(amblyObis$datasetName, "<br/>", amblyObis$catalogNumber, "<br/><a href=\"http://beta.iobis.org/explore/#/dataset/", amblyObis$resourceID, "\">OBIS dataset page</a>"), data = data.frame(lat = amblyObis$decimalLatitude, lng = amblyObis$decimalLongitude), radius = 3.5, weight = 0, fillColor = colors, fillOpacity = 1)

