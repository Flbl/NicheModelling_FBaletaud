#### Getting and formatting environmental variables ####

library(rmarkdown)
library(knitr)
library(dismo)
library(raster)
library(rgdal)
library(ncdf4)


## @knitr readingData

nc = nc_open(filename = "./Environment/global-analysis-forecast-phy-001-024.nc")
# print(nc)
# summary(nc)

ncRas = raster("./Environment/global-analysis-forecast-phy-001-024.nc", varname = "thetao")
ncRas
plot(ncRas, main = "Sea water potential temperature forecast (27-01-2017")

