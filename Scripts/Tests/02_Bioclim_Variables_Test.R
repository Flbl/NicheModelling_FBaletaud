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
nc = nc_open(filename = "./Environment/temp/rawData/3201152/CMEMS_GLO_001_024_TEMP_2007-01.nc")

nc = nc_open(filename = "./Environment/temp/rawData/3369745/CMEMS_GLO_001_024_TEMP_2016-03.nc")
nc = nc_open(filename = "./Environment/temp/rawData/5831837/CMEMS_GLO_001_024_TEMP_2011-11.nc")
nc

# Getting the temp data
ncTemps <- ncvar_get(nc,"thetao")
ncTemps

# Getting the long lat data
nclon <- ncvar_get(nc, "longitude")
nclat <- ncvar_get(nc, "latitude")


#Getting the time variable
ncTime <- ncvar_get(nc, "time")
ncTime
tunits <- ncatt_get(nc, "time", "units")
tunits
ncTimeDim <- dim(ncTime)
ncTimeDim


#Getting potential "NA" values in the data set
fillValues <- ncatt_get(nc, "thetao", "_FillValue")


#Splitting the time format for date conversion
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.numeric(unlist(tdstr)[2])
tday = as.numeric(unlist(tdstr)[3])
tyear = as.numeric(unlist(tdstr)[1])


# Converting the hours since 1950-01-01 to actual dates

#Conversion Nulle (se base sur du **days** since 1950-01-01 et pas moyen de changer dans les arguments)
t <- chron(ncTime, origin = c(day = tday, month = tmonth, year = tyear), format = "d/m/y")
t

#Conversion pas précise
as.Date(ncTime/24, origin="1950-01-01")

as.POSIXct(ncTime*60*60, origin="1950-01-01")

#Conversion précise
x <- as.POSIXct(ncTime, origin="1950-01-01", tz = "GMT") + as.difftime(ncTime,units="hours")
x[2500:3000]
#Autre conversion (avec package lubridate)
originDate <- ymd_hms("1950-01-01 00:00:00")
originDate + ncTime*3600





ncTemps[ncTemps == fillValues$value] <- NA
length(na.omit(ncTemps))




plot(earthGrid)
x <- ncTemps[,,1]
r <- raster(x, xmn = nclat[1], xmx = nclat[2], ymn = nclon[1], ymx = nclon[2], crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(r, add = TRUE)











##########################################
#using the second getCMEMS function


#getCMEMS
jan <- nc_open("./Scripts/getCMEMS/downs/global-analysis-forecast-phy-001-024_thetao_2013-01.nc")

# Read longitude & latitude
lon <- ncvar_get(jan, "longitude")
lat <- ncvar_get(jan, "latitude")

#Read the time
time_jan <- ncvar_get(jan, "time")

length(time_jan) == length(res)


#monthly

mon <- nc_open("downs/monthly_global-analysis-forecast-phy-001-024_thetao_2013-01.nc")

#Read the time
time_mon <- ncvar_get(mon, "time")

length(time_mon) == length(res_monthly)



