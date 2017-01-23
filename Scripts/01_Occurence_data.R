#################################################################
##                                                             ##
##            Getting and formatting occurence data            ##         
##                                                             ##
#################################################################

library(rgbif)
library(taxize)

# @knitr Download

ambly = name_lookup("Carcharhinus amblyrhynchos", rank = "species", return = "data")


name_suggest("Carcharhinus amblyrhynchos", rank = "species") #$key[1] for first suggested key (not necessarly the good one)
occ_search(scientificName = "Carcharhinus amblyrhynchos", limit = 20)
occ_count(taxonKey = 2418064, georeferenced = TRUE)

ambly = occ_search(scientificName = "Carcharhinus amblyrhynchos", fields = c("name", "basisOfRecord","protocol", return = "data"))
ambly = occ_search(scientificName = "Carcharhinus amblyrhynchos", return = "data")

gbifmap(ambly)
