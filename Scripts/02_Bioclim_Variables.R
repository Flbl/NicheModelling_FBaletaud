#################################################################
##                                                             ##
##      Getting and formatting the temperature variables       ##         
##                                                             ##
#################################################################


library(raster)
library(rgdal)
library(sp)
library(ncdf4)
library(ncdf.tools)
library(tools)
library(rgeos)
library(abind)

earthGrid <- raster("./data/interdata/earthGrid.tif")

eezNcGrid <- raster("./data/interdata/eezNcGrid.tif")


#########################################################################
# Code to get the temperature of cleandOccss cells only from copernicus #
#########################################################################

#1. Getting a list of unique cells 

getCellList <- function(cleandoccs){
  
  cellAll <-    unlist(lapply(cleandoccs,function(x){c(x$occs$cellNumber, x$abs$cellNumber)}))
  cellList <- unique(cellAll)
  
  return(cellList)
  
}

cellList <- getCellList(cleandOccs)


#### fixing for absences #########

# For the later use of the make file. The absence part needs to be fused with presence to make the script deal with both at the same time

cellListAbs1 <- unlist(lapply(cleandOccs,function(x){x$abs$CellNumber}))
cellListAbs1 <- unique(cellListAbs1)

x = unlist(lapply(cleandOccs,function(x){x$abs$decimalLongitude}))
y = unlist(lapply(cleandOccs,function(x){x$abs$decimalLatitude}))
df = data.frame(lon = x, lat = y)
df = unique(df[,c("lon","lat")])
dfmat <- split(as.matrix(df), seq(nrow(df)))
cellNumb <- sapply(dfmat, cellFromXY, object = earthGrid)
# test <- t(as.data.frame(sapply(cellNumb, getCellCentr)))
test <- t(sapply(cellNumb, getCellCentr))
testspt = SpatialPointsDataFrame(coords = test, data = as.data.frame(cellNumb), proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# writeOGR(testspt, dsn = "./data", layer = "testspt", driver = "ESRI Shapefile")
#########################
######## fixing absence data ########

cellListAbs <- cellNumb


#2. create a function returning cell centroid from cell number. (we know that those cells have occurrences)
getCellCentr <- function(cellID, raster = earthGrid){
  
  extOfCells <- coordinates(rasterFromCells(raster, cellID))
  
  extOfCells
  
}#eo getCellCentro

####################################################

#3. function that download temp data for a cell


clist <- as.list(cellList)
monthSeq <- as.list(append(paste0(0,seq(1:9)), as.character(10:12)))
yearSeq <- as.list(as.character(seq(0:9)+2000+6))

clistAbs <- as.list(cellListAbs)



## sourcing the function
source("./Scripts/getCMEMS/getCMEMS.R")

## parameters

### motu-client.py path
### this may chnage according to your computer configuration
motu_cl_lib <- "./Scripts/getCMEMS/libs/motu-client-python-master/src/python/motu-client.py"

### credentials (in cred.txt)
source("./Scripts/getCMEMS/cred.txt")



getCellTempData <- function(clist, monthSeq, yearSeq){
  
  outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMSabs/")
  
  lapply(clist, function(cellID){ # Inverser clist et yearSeq pour la prochaine fois
    
    lapply(yearSeq, function(year){
      
      lapply(monthSeq, function(month){
        
        prename = paste0("monthly_",cellID)
        
        #cell Extent
        cellExt <- getCellCentr(cellID)
        
        res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib ,
                                        log_cmems = log,
                                        pwd_cmems = pass,
                                        # Date 
                                        yyyystart=year,
                                        mmstart=month,
                                        # Area 
                                        xmin=as.character(cellExt[1,1]),
                                        xmax=as.character(cellExt[1,1]),
                                        ymin=as.character(cellExt[1,2]),
                                        ymax=as.character(cellExt[1,2]),
                                        zsmall="0.494", 
                                        zbig="1",
                                        #OutPath
                                        out_path = outDir,
                                        pre_name= prename)
        
        
        
      })
    })
  })
}


TempData <- getCellTempData(clist, monthSeq, yearSeq)
TempDataAbs <- getCellTempData(clistAbs, monthSeq, yearSeq)



#### Downloading missing data ####

fullCellFileListAbs <- as.vector(sapply(cellListAbs,function(cel){
  
  sapply(years,function(y){
    
    sapply(months,function(m){
      
      fName <- paste0("monthly_",cel,"global-analysis-forecast-phy-001-024_thetao_",y,"-",m,".nc")
      
    })
  })
}))

fullCellFileListAbs  

dlMissedTempData <- function(ncFile){
  
  out <- tryCatch({
    
    cat("#file :",ncFile,"\n")
    
    file <- nc_open(paste0("./data/rawdata/Environment/temp/CMEMS","/",ncFile))
    
    cat("Closing file")
    nc_close(file)
    
  }
  
  ,error = function(cond) {
    
    
    message("downloading data")
    
    param <- gsub("monthly_","", x = ncFile)
    param <- gsub("global-analysis-forecast-phy-001-024_thetao","",x = param)
    param <- gsub(".nc","", x = param)
    param <- unlist(strsplit(param, "[_-]+"))
    
    paramCell <- as.numeric(param[1])
    paramYear = param[2]
    paramMonth <- param[3]
    
    
    
    prename = paste0("monthly_",paramCell)
    
    #cell Extent
    cellExt <- getCellCentr(paramCell)
    
    res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib ,
                                    log_cmems = log,
                                    pwd_cmems = pass,
                                    # Date 
                                    yyyystart=paramYear,
                                    mmstart=paramMonth,
                                    # Area 
                                    xmin=as.character(cellExt[1,1]),
                                    xmax=as.character(cellExt[1,1]),
                                    ymin=as.character(cellExt[1,2]),
                                    ymax=as.character(cellExt[1,2]),
                                    zsmall="0.494", 
                                    zbig="1",
                                    #OutPath
                                    out_path = outDir,
                                    pre_name= prename)
    
  }
  
  )
  
  return(out)
  
}




TempDataMissed <- lapply(fullCellFileListAbs,dlMissedTempData)


cellListAll <- append(cellList, cellListAbs)











####################################################

# reading a nc files and storing them



getNCData <- function(ncFile){
  
  cat("#file :",ncFile,"\n")
  
  file <- nc_open(paste0("./data/rawdata/Environment/temp/CMEMS","/",ncFile))
  
  dat <- ncvar_get(file, varid = "thetao")
  
  res <- c(min(dat),max(dat),mean(dat))
  
  names(res) <- c("min","max","mean")
  
  nc_close(file)
  
  res
  
}#eo getNCData

# files <- list.files("./data/rawdata/Environment/temp/CMEMS")[1:12]

# resCell <- sapply(files,getNCData)



years <- as.character(2007:2016)
months <- c("01","02","03","04","05","06","07","08","09","10","11","12")


resCells <- lapply(cellListAll,function(cel){

  resYears <- lapply(years,function(y){

    sapply(months,function(m){

      fName <- paste0("monthly_",cel,"global-analysis-forecast-phy-001-024_thetao_",y,"-",m,".nc")

      getNCData(fName)

    })#eo lapply months

  })#eo lapply years

  resYears <- abind(resYears,along=3)

})#eo lapply cells

names(resCells) <- cellListAll





# 
# resCellsAbs <- lapply(cellListAbs,function(cel){
# 
#   resYears <- lapply(years,function(y){
# 
#     sapply(months,function(m){
# 
#       fName <- paste0("monthly_",cel,"global-analysis-forecast-phy-001-024_thetao_",y,"-",m,".nc")
# 
#       getNCData(fName)
# 
#     })#eo lapply months
# 
#   })#eo lapply years
# 
#   resYears <- abind(resYears,along=3)
# 
# })#eo lapply cells
# 
# names(resCellsAbs) <- cellListAbs

# resCellsAbs[[1]][1,,]





#Min
resCells[[1]][1,,]
resCells[[1]][,1,]
resCells[[1]][,,1:10]
mean(resCells[[1]][1,,])

# max

mean(resCells[[1]][2,,])

tempVar <- sapply(resCells, function(x){
  
   min = mean(x[1,,])
  
   max = mean(x[2,,])
   
   mean = mean(x[3,,])
   
   c(min,max,mean)
   # data.frame(min = min, max = max)
   
   
})



tempVar <- as.data.frame(t(tempVar))
colnames(tempVar) <- c("Tmin","Tmax","Tmean")
tempVar$cellNumber <- rownames(tempVar)
tempVar$Trange <- tempVar$Tmax - tempVar$Tmin
tempVar <- tempVar[,c(4,1,2,3,5)]
head(tempVar)
dim(tempVar)

write.csv(tempVar, "./data/calibdata/cellTemperatures.csv", row.names = FALSE)







######### Generating the predict data ##########


# Use maybe z-value from the raster layers to split data and make operations
# Then extract the temp value with the coordinates of the right raster (the one used in the whole study)
# NB : for the tool, maybe directly work from the extents from the first file downloaded (normally this one)


calTemp <- stack(rev(list.files("./data/rawdata/Environment/temp/New_Caledonia", full.names = TRUE))) # Rev because data is splitted in 2 files and we need the right date orders

dates <- seq(as.Date("2007-01-01"), as.Date("2016-12-31"), by="days")
dates <- as.character(dates)
names(calTemp) <- dates

# plot(calTemp[[3653]])
calTempdf <- as.data.frame(calTemp)

#get the date from the names of the layers and extract the year and month
indices <- unique(format(as.Date(names(calTemp), format = "X%Y.%m.%d"), format = "%Y.%m"))

# grep(indices[2], names(calTempdf))



monthlyMean <- sapply(indices, function(ind, data = calTempdf) {
  
  monthMean <- apply(data[,grep(ind, names(data))],1, mean)
  
  # monthMin <- apply(data[,grep(ind, names(calTempdf))],1, min)

  # monthMax <- apply(data[,grep(ind, names(calTempdf))],1, max)
  
  monthMean
}
)

monthlyMax <- sapply(indices, function(ind, data = calTempdf) {
  
  monthMean <- apply(data[,grep(ind, names(data))],1, max)
  
  # monthMin <- apply(data[,grep(ind, names(calTempdf))],1, min)
  
  # monthMax <- apply(data[,grep(ind, names(calTempdf))],1, max)
  
  monthMean
}
)


monthlyMin <- sapply(indices, function(ind, data = calTempdf) {
  
  monthMean <- apply(data[,grep(ind, names(data))],1, min)
  
  # monthMin <- apply(data[,grep(ind, names(calTempdf))],1, min)
  
  # monthMax <- apply(data[,grep(ind, names(calTempdf))],1, max)
  
  monthMean
}
)

meanOfMax <- apply(monthlyMax, 1, mean)
meanOfMin <- apply(monthlyMin, 1, mean)
meanOfMean <- apply(monthlyMean, 1, mean)

ncTempVar <- data.frame(Tmin = meanOfMin, Tmax = meanOfMax, Tmean = meanOfMean)
ncTempVar$Trange <- ncTempVar$Tmax - ncTempVar$Tmin
ncTempVar <- cbind(as.data.frame(coordinates(calTemp)), ncTempVar)
ncTempVar <- na.omit(ncTempVar)

ncTempVarSpt <- SpatialPointsDataFrame(
  coords = cbind(ncTempVar$x, ncTempVar$y),
  data = ncTempVar,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
)


## Extracting the values to the true centroids of cell of grid used for the rest of the study

trueNcTempVar <- data.frame(datToRemove = getValues(eezNcGrid))

#True min
TminRas <- raster(extent(calTemp))
res(TminRas) <- res(calTemp)
projection(TminRas) <- proj4string(calTemp)
origin(TminRas) <- origin(calTemp)
TminRas <- rasterize(ncTempVarSpt,field = "Tmin", TminRas)
# plot(TminRas)
trueNcTempVar$Tmin <- extract(TminRas, coordinates(eezNcGrid))

#True max
TminRas <- raster(extent(calTemp))
res(TminRas) <- res(calTemp)
projection(TminRas) <- proj4string(calTemp)
origin(TminRas) <- origin(calTemp)
TminRas <- rasterize(ncTempVarSpt,field = "Tmax", TminRas)
# plot(TminRas)
trueNcTempVar$Tmax <- extract(TminRas, coordinates(eezNcGrid))


#True mean
TminRas <- raster(extent(calTemp))
res(TminRas) <- res(calTemp)
projection(TminRas) <- proj4string(calTemp)
origin(TminRas) <- origin(calTemp)
TminRas <- rasterize(ncTempVarSpt,field = "Tmean", TminRas)
# plot(TminRas)
trueNcTempVar$Tmean <- extract(TminRas, coordinates(eezNcGrid))


#True range
TminRas <- raster(extent(calTemp))
res(TminRas) <- res(calTemp)
projection(TminRas) <- proj4string(calTemp)
origin(TminRas) <- origin(calTemp)
TminRas <- rasterize(ncTempVarSpt,field = "Trange", TminRas)
# plot(TminRas)
trueNcTempVar$Trange <- extract(TminRas, coordinates(eezNcGrid))



# test[which(is.na(test$datToRemove) == TRUE),] <- NA

#####
trueNcTempVar[which(is.na(trueNcTempVar$datToRemove) == TRUE),] <- NA # Removing data outside EEZ
trueNcTempVar <- cbind(coordinates(eezNcGrid), trueNcTempVar)
trueNcTempVar <- na.omit(trueNcTempVar)
trueNcTempVar$datToRemove <- NULL



# @knitr WriteToRast

writeVarToRast <- function(s, map = eezNcGrid){
  
  sNoCoords <-  subset(s, select = -c(x,y))
  
  substrateSpatial <-  SpatialPointsDataFrame(
    coords = cbind(s$x, s$y),
    data = sNoCoords,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  featRast <- raster(extent(map))
  res(featRast) <- res(map)
  projection(featRast) <- proj4string(map)
  origin(featRast) <- origin(map)
  
  lapply(names(sNoCoords), function(x){
    
    rast <- featRast
    
    rast <- rasterize(substrateSpatial ,field = x, featRast)
    
    writeRaster(rast, filename = paste0("./data/predictdata/BEM", x, ".tif"), overwrite = TRUE)
    
    
  })
  
  
  
}# eo writeVarToRast


writeVarToRast(trueNcTempVar)


# plot(raster("./data/predictdata/BEM/Tmean.tif"))


