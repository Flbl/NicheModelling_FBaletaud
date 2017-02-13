##########################################
## Formatting the temperature variables ##
##########################################



library(raster)
library(rgdal)
library(sp)
library(ncdf4)
library(tools)
library(rgeos)
library(chron)

# reading a cell nc file

nc = nc_open(filename = "./Environment/temp/rawData/4622441/CMEMS_GLO_001_024_TEMP_2007-01.nc")
ncTemps <- ncvar_get(nc,"thetao")

nclon <- ncvar_get(nc, "longitude")
nclat <- ncvar_get(nc, "latitude")

ncTime <- ncvar_get(nc, "time")
tunits <- ncatt_get(nc, "time", "units")
ncTimeDim <- dim(ncTime)
fillValues <- ncatt_get(nc, "thetao", "_FillValue")


tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.numeric(unlist(tdstr)[2])
tday = as.numeric(unlist(tdstr)[3])
tyear = as.numeric(unlist(tdstr)[1])
t <- chron(ncTime, origin = c(day = tday, month = tmonth, year = tyear), format = "d/m/y")
t[1000:1200]

ncTemps[ncTemps == fillValues$value] <- NA
length(na.omit(ncTemps))

nc = nc_open(filename = "./Environment/temp/rawData/4622441/CMEMS_GLO_001_024_TEMP_2007-01.nc")




