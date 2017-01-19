#### Getting and formatting environmental variables ####

library(rmarkdown)
library(knitr)
library(dismo)
library(raster)
library(rgdal)
library(sp)
library(ncdf4)
library(tools)

## @knitr temperature

nc = nc_open(filename = "./Environment/global-analysis-forecast-phy-001-024.nc")
# print(nc)
# summary(nc)

ncRas = raster("./Environment/global-analysis-forecast-phy-001-024.nc", varname = "thetao")
ncRas
plot(ncRas, main = "Sea water potential temperature forecast (27-01-2017")



## @knitr seafloor

shps = dir("./Environment/global-seafloor-geomorphic-features-map", "*.shp")
shps = file_path_sans_ext(shps)


for(shp in shps) {
  assign(shp, readOGR("./Environment/global-seafloor-geomorphic-features-map",layer = shp))
}

# shpv = as.vector(shps)
# 
# for(shpras in shp){
#   shpras = raster(shp)
#   res(shpras) = 1/12
#   shpras = rasterize(shp, field = 1, shapras)
# }

Seamountsras = raster(Seamounts)
res(Seamounts) = 1/12
Seamountsras = rasterize(Seamounts, field = "Geomorphic", Seamountsras)

plot(Seamounts, add = TRUE)
plot(Rift_valleys, add = TRUE)         
plot(Plateaus, add = TRUE)
