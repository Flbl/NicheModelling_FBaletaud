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

getCellTempData <- function(cellID){
  
  
  #cell Extent
  cellExt <- getCellCentr(cellID)

  ## sourcing the function
  source("./Scripts/getCMEMS/getCMEMS.R")

  ## parameters

  ### motu-client.py path
  ### this may chnage according to your computer configuration
  motu_cl_lib <- "./Scripts/getCMEMS/libs/motu-client-python-master/src/python/motu-client.py"

  ### output dir
  outDir <- paste0("./data/rawdata/Environment/temp/CMEMS",cellID,"/")

  ### credentials (in cred.txt)
  source("./Scripts/getCMEMS/cred.txt")




res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib , 
                                out_path = outDir,
                                log_cmems = log,
                                pwd_cmems = pass,
                                # Date 
                                yyyystart="2013",
                                mmstart="01",
                                # Area 
                                xmin=as.character(cellExt[1,1]),
                                xmax=as.character(cellExt[1,1]),
                                ymin=as.character(cellExt[1,2]),
                                ymax=as.character(cellExt[1,2]))
}




























