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
  
 lapply(yearSeq, function(year){
  
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





































