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
  
  cellAll <-    unlist(lapply(cleandOccs,function(x){c(x$occs$cellNumber, x$abs$cellNumber)}))
  cellList <- unique(cellAll)
  
  return(cellList)
  
}

cellList <- getCellList(cleandOccs)



#2. create a function returning cell centroid from cell number. (we know that those cells have occurrences)
getCellCentr <- function(cellID, raster = earthGrid){
  
  extOfCells <- coordinates(rasterFromCells(raster, cellID))
  
  extOfCells
  
}#eo getCellCentro

####################################################

#3. function that download temp data for a cell


clist <- as.list(cellList)
monthSeq <- append(paste0(0,seq(1:9)), as.character(10:12))
monthSeq <- as.list(monthSeq)
yearSeq <- as.list(as.character(seq(0:9)+2000+6))




## sourcing the function
source("./Scripts/getCMEMS/getCMEMS.R")

## parameters

### motu-client.py path
### this may chnage according to your computer configuration
motu_cl_lib <- "./Scripts/getCMEMS/libs/motu-client-python-master/src/python/motu-client.py"

### credentials (in cred.txt)
source("./Scripts/getCMEMS/cred.txt")



getCellTempData <- function(clist, monthSeq, yearSeq){
  
  outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/")
  
  lapply(yearSeq, function(year){ # Inverser clist et yearSeq pour la prochaine fois
    
    lapply(monthSeq, function(month){
      
      lapply(clist, function(cellID){
        
        prename= paste0("monthly_",cellID)
        
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


resCells <- lapply(cellList,function(cel){
  
  resYears <- lapply(years,function(y){
    
    sapply(months,function(m){
      
      fName <- paste0("monthly_",cel,"global-analysis-forecast-phy-001-024_thetao_",y,"-",m,".nc")
      
      getNCData(fName)
      
    })#eo lapply months
    
  })#eo lapply years
  
  resYears <- abind(resYears,along=3)
  
})#eo lapply cells

names(resCells) <- cellList

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
   
   c(min,max)
   # data.frame(min = min, max = max)
   
   
})

# class(as.data.frame(tempVar))
# 
# df <- as.data.frame(tempVar)
