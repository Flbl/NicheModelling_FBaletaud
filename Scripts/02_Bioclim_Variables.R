#################################################################
##                                                             ##
##      Getting and formatting the temperature variables       ##         
##                                                             ##
#################################################################


library(raster)
library(rgdal)
library(sp)
library(ncdf4)
library(tools)
library(rgeos)

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



# reading cell files and storing them

cellFileList <- dir("./data/rawdata/Environment/temp/CMEMS")

#Function (to do) generating variable for one cell

getCellTempVar <- function(clist, cellFileList){
  
  
  cellData <- lapply(clist, function(cellID, cf = cellFileList){
  
  
  fname <- cf[grep(as.character(cellID),cf)]
  
  fname <- paste0("./data/rawdata/Environment/temp/CMEMS","/", fname)
  
  cellFiles <- lapply(fname, nc_open)
  
  celldata <- lapply(cellFiles, ncvar_get, varid = "thetao")
  
  
  #annual mean
  meanTemp <-  mean(unlist(celldata))
  
  #minimum
  minTemp <- min(unlist(celldata))
  minTemp <- lapply(celldata, mean)
  
  #Max
  maxTemp <- max(unlist(celldata))
  
  #Annual Range
  tempRange <- maxTemp - meanTemp
  
  
  df <- data.frame(cell = cellID, MAT = meanTemp, MINT = minTemp, MAXT = maxTemp, TRAN = tempRange)
  
  # rownames(df) <- cellID
  
  df
  
  
  }
)

  
  
}























