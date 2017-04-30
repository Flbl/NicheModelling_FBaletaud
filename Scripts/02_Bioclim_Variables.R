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
tempVar$CellNumber <- rownames(tempVar)
tempVar$Trange <- tempVar$Tmax - tempVar$Tmin
tempVar <- tempVar[,c(4,1,2,3,5)]
head(tempVar)
dim(tempVar)

write.csv(tempVar, "./data/calibdata/cellTemperatures")
