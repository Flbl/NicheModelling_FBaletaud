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

nc = nc_open(filename = "./Environment/METOFFICE-GLO-SST-L4-NRT-OBS-SST-MON-V2_1484923518762.nc")
print(nc)
# summary(nc)

ncRas = raster("./Environment/METOFFICE-GLO-SST-L4-NRT-OBS-SST-MON-V2_1484923518762.nc", varname = "analysed_sst")
ncRas
plot(ncRas, main = "Analysed Sea Surface Temperature")



## @knitr seafloor

shps = dir("./Environment/global-seafloor-geomorphic-features-map", "*.shp")
shps = file_path_sans_ext(shps)


for(shp in shps) {
  assign(shp, readOGR("./Environment/global-seafloor-geomorphic-features-map",layer = shp))
}

CoralShp = readOGR("./Environment/DataPack-14_001_WCMC008_CoralReef2010_v1_3/01_Data", layer = "14_001_WCMC008_CoralReef2010_v1_3")

# shpv = as.vector(shps)
# 
# for(shpras in shp){
#   shpras = raster(shp)
#   res(shpras) = 1/12
#   shpras = rasterize(shp, field = 1, shapras)
# }



Seamountsras = raster(Seamounts)
res(Seamountsras) = 1/12
Seamountsras = rasterize(Seamounts, field = "Geomorphic", Seamountsras)

plot(Seamounts, add = TRUE)
plot(Rift_valleys, add = TRUE)         
plot(Plateaus, add = TRUE)
