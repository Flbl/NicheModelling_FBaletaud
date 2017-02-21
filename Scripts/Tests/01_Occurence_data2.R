#################################################################
##                                                             ##
##            Getting and formatting occurence data            ##         
##                                                             ##
#################################################################

library(rgbif)
library(taxize)
library(taxizesoap)
library(devtools)
library(robis)
library(leaflet)
library(tibble)
library(rgdal)
library(dplyr)
library(sp)

library(ggmap)
library(dplyr)
library(ncdf4)
library(raster)
library(stringr)
library(readr)
library(pryr)
library(marmap)
library(lubridate)
library(broom)
library(mregions)
library(rgeos)



# Reading/creating the earth GRID a transformer en fonction

earthGrid<-nc_open("./Environment/temp/earthgrid/global-analysis-forecast-phy-001-024.nc") # Open any nc file from the bioclimatic data used


earthGrid<-ncvar_get(earthGrid,"thetao")



earthGrid[(is.na(earthGrid))]<-1000
earthGrid[(earthGrid<1000)]<-NA
earthGrid[(earthGrid == 1000)]<-1
earthGrid<-earthGrid[,order(c(1:ncol(earthGrid)),decreasing=T)] # Inverting coordinates to replace them right

earthGrid <- raster(nrows = 2041, ncols = 4320, xmn = -180, xmx = 179.9167, ymn = -80, ymx = 90, crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0", vals = as.vector(earthGrid)) # Check the coordinates of the frame on copernicus

earthRegion <- readOGR(dsn = "./Environment/GSHHS_region", layer = "GSHHS_f_L1")







getCleanedOcc = function(my_sp, qc = c(1:7,10:19,21:30)){
  
  # CheckSp function
  checkSp <- gsub("_", " ",my_sp)
  
  
  checkSp_aphia <- get_wormsid(searchterm = checkSp, accepted = FALSE)
  checkSp_taxo <- worms_records(ids = checkSp_aphia, marine_only = TRUE)
  
  if(is.na(checkSp_aphia)){
    stop("Species not found")
  }
  
  if(checkSp_taxo$status == "unaccepted"){
    print("invalid name : assigning latest name")
    my_sp <- checkSp_taxo$valid_name
  }
  if(checkSp_taxo$status == "accepted"){
    print("Valid name")
    my_sp <- checkSp
    
  }
  review = data.frame(Entered_Name = checkSp, Reason = checkSp_taxo$unacceptreason, Accepted_Name = my_sp)
  print(review)
  
  
  
  
  # OBISSP function
  
  print("total number of records available")
  my_sp_summ = checklist(my_sp)
  print(my_sp_summ$records)
  
  print("Getting data through Quality Flags :")
  
  OBIS_occs = occurrence(my_sp, qc = qc)
  
  print("Keeping only occurrences after 1980")
  OBIS_occs = OBIS_occs[OBIS_occs$yearcollected >= 1980,]
  
  # if (nrow(OBIS_occs) <= 2) stop("Not enough records")
  
  #assign(paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2],"OBIS_occs", sep = "", collapse = ""), OBIS_occs, .GlobalEnv)
  
  
  
  
  
  # Merging the occurrences from OBIS and New Caledonian dataset
  
  ## Adding occurrence column to OBIS dataset
  if (nrow(OBIS_occs) >= 1) OBIS_occs$occurrence = 1
  
  ## Reading/formatting the dataset
  NCRecords <- read.table("./Biodiversity/dataset_fitted.csv", sep = ";", dec = ",", header = TRUE) 
  print("Creating/Formatting the NC dataset")
  match = paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2], sep = "", collapse = "")
  matchCol = grep(match,names(NCRecords))
  names(NCRecords)[matchCol] <- my_sp
  
  NC_occs <- NCRecords[,c(2,3, matchCol)]
  NC_occs <- NC_occs[NC_occs[,3] >= "1",]
  colnames(NC_occs) <- c("decimalLatitude","decimalLongitude","individualCount")
  NC_occs$occurrence = 1
  NC_occs$datasetName <- "APEX"
  rownames(NC_occs) <- NULL # Reassign the rownames of NC_occs
  
  
  #assign(paste(substring(my_sp, 1,1), "_" ,unlist(strsplit(my_sp, " "))[2],"NC_occs", sep = "", collapse = ""), NC_occs, .GlobalEnv)
  
  
  ## merge OBIS and NC
  
  print("Merging OBIS and NC datasets")
  OccsRaw <- merge(OBIS_occs, NC_occs, all = TRUE)

  Occs <- OccsRaw[, names(NC_occs)]

  
  
  
  
  
  
  # Manual check of the dataset
  
  print("deleting coordinates duplicates")
  
  if(any(is.na(Occs$decimalLatitude) & is.na(Occs$decimalLongitude))) Occs = Occs[-which(is.na(Occs$decimalLatitude) & is.na(Occs$decimalLongitude)),]
  
  Occs = distinct(Occs, decimalLongitude, decimalLatitude, .keep_all = TRUE)
  
  print("replacing NA with 1 to rows without individual count (fisheries)")
  
  Occs[is.na(Occs$individualCount),3] <- 1
  
  rownames(Occs) <- NULL # Reassign the rownames of NC_occs
  
  
  
  
  
  
  
  
  # checking for points on earth using the bioclimatic grid
  
  print("Check points on earth using bioclimatic grid")
  
  
  ## creating spatial points for our occurrences
  
  occsSpatial = SpatialPointsDataFrame(
    coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
    data = Occs,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  )
  
  # plot(occsSpatial, add = TRUE)
  
  Occs <- Occs[-which(raster::extract(earthGrid,occsSpatial) == 1),]
  
  # writeOGR(occsSpatial, dsn = "./environment/temp", layer = "OccsSpatialNferru", driver = "ESRI Shapefile")
  
  
  rownames(Occs) <- NULL # Reassign the rownames of Occs
  
  
  
  
  
  
  # Double checking with the GSHHS "special islets" regions shapefile
  
  print("Check points on earth using GSHHS shapefile")
  
  
  occsSpatial = SpatialPointsDataFrame(
    coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
    data = Occs,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  )
  
  # earthRegion <- readOGR(dsn = "./Environment/GSHHS_region", layer = "GSHHS_f_L1")
  
  is_occ_terrest <- gIntersects(occsSpatial,earthRegion, byid = TRUE, returnDense = FALSE)
  
  occ_terrest <- as.numeric(names(unlist(is_occ_terrest)))
  
  Occs <- Occs[-occ_terrest,]
  
  rownames(Occs) <- NULL
  
  # Adding the species name
  Occs$species <- gsub(" ", "_",my_sp)
  
  #Getting the nrow of unique cell number/cell with occurrence
  # = Remove all but one occurrence per cell = convert occurrences to cells at 1/12
  cat("Getting the cell numbers of occurrence\n")
  
  occsSpatial = SpatialPointsDataFrame(
    coords = cbind(Occs$decimalLongitude, Occs$decimalLatitude),
    data = Occs,
    proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  
  
  OccsRast <- raster(extent(earthGrid))
  
  res(OccsRast) <- res(earthGrid)
  
  projection(OccsRast) <- proj4string(occsSpatial)
  
  origin(OccsRast) <- origin(earthGrid)
  
  OccsRast <- rasterize(occsSpatial,field = "occurrence", OccsRast)
  
  coords <- cbind(Occs$decimalLongitude, Occs$decimalLatitude)
  
  Occs$cellNumber <- cellFromXY(OccsRast, coords)
  
  
  
  
  return(Occs)
  
  
}#eo GetCleanedOcc



generateAbs <- function(Occs){
  
  # Creating pseudo absences using IUCN range map outer limits
  print("Getting IUCN range map data")
  
  ## Generate auto paths to species shapefile
  my_sp <- Occs$species[1]
  
  pathSp = paste("./Biodiversity/iucn/", my_sp, sep = "", collapse = "")
  
  layerSp = sub(".shp", "", list.files(path = pathSp, full.names = FALSE, pattern = ".shp")[1])
  
  
  ## Reading the shapefile of the species
  cat("reading the IUCN shapefile\n")
  
  iucn = readOGR(dsn = pathSp, layer = layerSp)
  
  
  ## Converting to raster
  cat("Converting shapefile to raster\n")
  
  rasIucnOccs <- raster(extent(iucn))
  
  res(rasIucnOccs) <- res(earthGrid)
  
  projection(rasIucnOccs) <- proj4string(iucn)
  
  origin(rasIucnOccs) <- origin(earthGrid)
  
  rasIucnOccs <- rasterize(iucn, field = "PRESENCE", rasIucnOccs) # assigning 1 for cells within the range
  
  cat("merging earthGrid and IUCN rangemap raster\n")
  
  antiAbs <- merge(earthGrid, rasIucnOccs)
  
  # antiAbs[is.na(antiAbs)] <- 2
  
  antiAbs <- reclassify(antiAbs, cbind(NA,2)) # Assigning 2 for cells outside the range
  
  
  
  
  
  
  
  
  # Sampling the cells with the length of OccCells in the cells exterior to IUCN range and earthGrid
  
  OccsCells <- unique(Occs$cellNumber)
  
  cat("sampling cells outside species range and earth\n")
  
  sampleClasses <- function(r, n)  {
    
      cellVal <- which(t(as.matrix(r)) == 2) # get All cells for class 2 (=Abs)
      
      samples <- sample(cellVal, n) # sample class's cell number
      
      return(samples)
    }
    
    
  pseudoAbsCells <- sampleClasses(antiAbs, length(OccsCells))
  
  pseudoAbs <- as.data.frame(xyFromCell(antiAbs, pseudoAbsCells))
  
  colnames(pseudoAbs) <- c("decimalLongitude","decimalLatitude")
  
  pseudoAbs$CellNumber <- pseudoAbsCells
  
  pseudoAbs$species <- Occs$species[1]
  
  pseudoAbs

  
} # eo generateAbs




getOccsAbs <- function(my_sp){
  
  cat("entering getcleanedOcc\n")
  Occs <- getCleanedOcc(my_sp)
  
  cat("entering getOccsAbs\n")
  Abs <-generateAbs(Occs)
  
  list(occs=Occs, abs=Abs)
  
  
}




species <- c("Carcharhinus_amblyrhynchos","Carcharhinus_melanopterus","Triaenodon_obesus")


cleandOccs <- lapply(species,getOccsAbs)

names(cleandOccs) <- species


# cleandOccs[["Carcharhinus_amblyrhynchos"]] # <=> cleandOccs$Carcharhinus_amblyrhynchos
# 
# cleandOccsMerged <- do.call(rbind,cleandOccs) # -> rbind tous les tableaux de la liste







####################################
# recup env pour les cellules eelctionnées
####################################

#1. obtenir une liste de cellules uniques

getCellList <- function(cleandoccs){

  cellAll <-    unlist(lapply(cleandOccs,function(x){c(x$occs$cellNumber, x$abs$cellNumber)}))
  cellList <- unique(cellAll)

  return(cellList)

}

cellList <- getCellList(cleandOccs)

# #2. create a function returning cell extent from cell number
# getCellExtent <- function(cellList, raster = earthGrid){
# 
#   extOfCells <- lapply(lapply( X = cellList, rasterFromCells, x = raster), extent)
# 
#   names(extOfCells) <- cellList
# 
#   extOfCells
# 
# }#eo getCellExtent


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
getCellTempData(cellList[1])


lapply(cellList,getCellTempData)























getCellTempData <- function(cellID){
  
  #out dir
  patName <- paste0("/home/florian/Travail_FBaletaud/Environment/temp/rawData/",cellID)
  dir.create(patName)
  outDir <- paste0("D:\\Florian\\Documents\\Universite\\M1SBM\\Stage\\Travail_FBaletaud\\Environment\\temp\\rawData\\",cellID,"\\")
  
  #cell Extent
  cellExt <- getCellExtent(cellID)
  
  
  # cellID <- cellExt[[1]]
  
  #call python stuff
  
  ### python.exe path
  myPythonPath <- "D:\Florian\Documents\Universite\M1SBM\Stage\Python27\python.exe"
  
  ## sourcing the function
  source("D:/Florian/Documents/Universite/M1SBM/Stage/CMEMS3567_GLO_Daily_by_Month_R/getCMEMS.R")
  
  ## parameters
  
  ### motu-client.py path
  motu_cl_lib <- "D:/Florian/Documents/Universite/M1SBM/Stage/Travail_FBaletaud/Scripts/CMEMS3567_GLO_Daily_by_Month_R/libs/motu-client-python-master/src/python/motu-client.py"
  
  
  ### call the function
  getCMEMS(scriptPath="D:/Florian/Documents/Universite/M1SBM/Stage/Travail_FBaletaud/Scripts/CMEMS3567_GLO_Daily_by_Month_R/libs/CMEMS3567_GLO_Daily_by_Month_CallFromR.py",
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







#####################################


# Cambly <- getOccsAbs("Carcharhinus_amblyrhynchos")
# 
# CamblyAbs <- Cambly[[2]]



# 
# getCellsExtent <- function(occs, abs){
# 
# 
#   extOfCells <- lapply(lapply( X = OccsCells, rasterFromCells, x = OccsRast), extent)
# 
#   # extOfCells[[1]]
# 
#   extOfCells
# 
# }
# 
# 
# Camblycells_ext <- getCellsExtent(Occs)




# Plotting for verification
leaflet() %>%
  addProviderTiles("Esri.WorldImagery") %>%
  addCircleMarkers(data = data.frame(lat = Occs$decimalLatitude[21], lng = Occs$decimalLongitude[21]), radius = 3.5, weight = 0, fillColor = "orange", fillOpacity = 1)




