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

nc = nc_open(filename = "./Environment/temp/A20160012016366.L3m_YR_SST4_sst4_4km.nc")
print(nc)
# summary(nc)

ncRas = raster("./Environment/temp/A20160012016366.L3m_YR_SST4_sst4_4km.nc", varname = "sst4")
ncRas
plot(ncRas, main = "Sea Surface Temperature")

## Temperature NC

TNC <- read.table("./environment/temp/Temp_NC_clean.csv" , sep = ";" , dec = "," , header = TRUE)
head(TNC)
TNC = TNC[,c(2,3,1)]
TNCRast <- rasterFromXYZ(TNC)
crs(TNCRast) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(TNCRast)




TNCSpatial <- SpatialPointsDataFrame(
  coords = cbind(TNC$Longitude, TNC$Latitude),
  data = TNC,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)

TNCSpatial
TNCGrid <- TNCSpatial
gridded(TNCGrid) <- TRUE
TNCGrid <- as(TNCGrid, "SpatialGridDataFrame")
TNCPix <- as(TNCGrid,"SpatialPixelsDataFrame")
TNCRast <- raster(TNCPix)
TNCras <- raster(TNCGrid)
TNCras <- rasterize(TNCSpatial, TNCras)



## @knitr Terrest

Region <-readOGR("./Environment/GSHHS_region","GSHHS_f_L1")

RegionDissolved <- gUnionCascaded(Region)

d <- data.frame(c(1))

ecorDissolved <- SpatialPolygonsDataFrame(ecorDissolved,data=d)

writeOGR(ecorDissolved, ".", "GSHHS_f_L1_dissolved", driver="ESRI Shapefile")


## @knitr seafloor

shps = dir("./Environment/global-seafloor-geomorphic-features-map", "*.shp")
shps = file_path_sans_ext(shps)


for(shp in shps) {
  assign(shp, readOGR("./Environment/global-seafloor-geomorphic-features-map",layer = shp))
}

CoralShp = readOGR("./Environment/WCMC008_CoralReef2010_v1_3/01_Data", layer = "14_001_WCMC008_CoralReef2010_v1_3")

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


