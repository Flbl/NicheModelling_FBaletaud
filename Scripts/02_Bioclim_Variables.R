##########################################
## Formatting the temperature variables ##
##########################################



library(raster)
library(rgdal)
library(sp)
library(ncdf4)
library(tools)
library(rgeos)


# reading a 
nc = nc_open(filename = "./Environment/temp/A20160012016366.L3m_YR_SST4_sst4_4km.nc")
