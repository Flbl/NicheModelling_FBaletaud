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



#2. create a function returning cell extent from cell number
getCellExtent <- function(cellID, raster = earthGrid){
  
  extOfCells <- extent(rasterFromCells(raster, cellID))
  
  extOfCells
  
}#eo getCellExtent



#3. function that download temp data for a cell


getCellTempData <- function(cellID){
  
  #out dir
  patName <- paste0("./Environment/temp/rawData/",cellID)
  dir.create(patName)
  outDir <- paste0("./Environment/temp/rawData/",cellID,"/")
  
  #cell Extent
  cellExt <- getCellExtent(cellID)
  
  # cellID <- cellExt[[1]]
  
  
  
  #call python stuff
  
  ## sourcing the function
  source("./Scripts/CMEMS3567_GLO_Daily_by_Month_R/getCMEMS.R")
  
  ## parameters
  
  myPythonPath <- "python"
  
  ## motu-client.py path
  ### this may change according to your computer configuration
  motu_cl_lib <- "./Scripts/CMEMS3567_GLO_Daily_by_Month_R/libs/motu-client-python-master/src/python/motu-client.py"
  
  
  ### call the function
  getCMEMS(scriptPath="./Scripts/CMEMS3567_GLO_Daily_by_Month_R/libs/CMEMS3567_GLO_Daily_by_Month_CallFromR.py",
           python=myPythonPath,
           motu_cl = motu_cl_lib,
           log_cmems="fbaletaud",
           pwd_cmems="FlorianCMEMS2017",
           #Date
           yyyystart="2006",
           mmstart="12",
           yyyyend="2016",
           mmend="12",
           hh=" 12:00:00",
           dd="31",
           # Area 
           xmin=as.character(cellExt@xmin),
           xmax=as.character(cellExt@xmax),
           ymin=as.character(cellExt@ymin),
           ymax=as.character(cellExt@ymax),
           zmin="0.49", 
           zmax="0.50", 
           # Variables 
           table_var_cmd = "thetao",
           table_data_type = "TEMP_",
           # Output files 
           out_path =  outDir, #Make sure to end your path with "/" 
           pre_name= "CMEMS_GLO_001_024_")
  
  
}#eo getCellTempData



#test
# getCellTempData(cellList[1])

#get Data
lapply(cellList,getCellTempData)








# reading a 
nc = nc_open(filename = "./Environment/temp/A20160012016366.L3m_YR_SST4_sst4_4km.nc")











