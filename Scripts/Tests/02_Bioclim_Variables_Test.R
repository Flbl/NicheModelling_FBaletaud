##########################################
## Formatting the temperature variables ##
##########################################



library(raster)
library(rgdal)
library(sp)
library(ncdf4)
library(tools)
library(rgeos)
library(chron)

# reading a cell nc file

nc = nc_open(filename = "./data/rawdata/Environment/temp/CMEMS/monthly_3734920global-analysis-forecast-phy-001-024_thetao_2012-11.nc")
nc = nc_open(filename = "./Environment/temp/rawData/3201152/CMEMS_GLO_001_024_TEMP_2007-01.nc")

nc = nc_open(filename = "./Environment/temp/rawData/3369745/CMEMS_GLO_001_024_TEMP_2016-03.nc")
nc = nc_open(filename = "./Environment/temp/rawData/5831837/CMEMS_GLO_001_024_TEMP_2011-11.nc")
nc

# Getting the temp data
ncTemps <- ncvar_get(nc,"thetao")
ncTemps

# Getting the long lat data
nclon <- ncvar_get(nc, "longitude")
nclat <- ncvar_get(nc, "latitude")


#Getting the time variable
ncTime <- ncvar_get(nc, "time")
ncTime
tunits <- ncatt_get(nc, "time", "units")
tunits
ncTimeDim <- dim(ncTime)
ncTimeDim


#Getting potential "NA" values in the data set
fillValues <- ncatt_get(nc, "thetao", "_FillValue")


#Splitting the time format for date conversion
tustr <- strsplit(tunits$value, " ")
tdstr <- strsplit(unlist(tustr)[3], "-")
tmonth = as.numeric(unlist(tdstr)[2])
tday = as.numeric(unlist(tdstr)[3])
tyear = as.numeric(unlist(tdstr)[1])


# Converting the hours since 1950-01-01 to actual dates

#Conversion Nulle (se base sur du **days** since 1950-01-01 et pas moyen de changer dans les arguments)
t <- chron(ncTime, origin = c(day = tday, month = tmonth, year = tyear), format = "d/m/y")
t

#Conversion pas précise
as.Date(ncTime/24, origin="1950-01-01")

as.POSIXct(ncTime*60*60, origin="1950-01-01")

#Conversion précise
x <- as.POSIXct(ncTime, origin="1950-01-01", tz = "GMT") + as.difftime(ncTime,units="hours")
x[2500:3000]
#Autre conversion (avec package lubridate)
originDate <- ymd_hms("1950-01-01 00:00:00")
originDate + ncTime*3600





ncTemps[ncTemps == fillValues$value] <- NA
length(na.omit(ncTemps))




plot(earthGrid)
x <- ncTemps[,,1]
r <- raster(x, xmn = nclat[1], xmx = nclat[2], ymn = nclon[1], ymx = nclon[2], crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(r, add = TRUE)











##########################################
#using the second getCMEMS function


#getCMEMS
jan <- nc_open("./Scripts/getCMEMS/downs/global-analysis-forecast-phy-001-024_thetao_2013-01.nc")

# Read longitude & latitude
lon <- ncvar_get(jan, "longitude")
lat <- ncvar_get(jan, "latitude")

#Read the time
time_jan <- ncvar_get(jan, "time")

length(time_jan) == length(res)


#monthly

mon <- nc_open("./Scripts/getCMEMS/downs/monthly_global-analysis-forecast-phy-001-024_thetao_2013-01.nc")

# Read values
temp_mon <- ncvar_get(jan,"thetao")

lon <- ncvar_get(mon, "longitude")
lat <- ncvar_get(mon, "latitude")

#Read the time
time_mon <- ncvar_get(mon, "time")

length(time_mon) == length(res_monthly)


# Convert to raster

# monrast <- raster(temp_mon[,,1], xmn = lon[1], xmx = lon[2], ymn = lat[1], ymx = lat[2], crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
xyz <- data.frame(x=lon, y = lat, z = temp_mon[1])
monrast <- rasterFromXYZ(xyz, res = res(earthGrid), crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
plot(monrast)

#occ point from this data
occPt <- SpatialPoints(coords = cbind(cleandOccs[[1]][[1]][1,2], cleandOccs[[1]][[1]][1,1]),
              proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(occPt, add = TRUE)







## Trying to make getCellTempsData to work

getCellTempData <- function(cellID,
                            python = "python",  #path to python executable
                            motu_cl = "libs/motu-client-python-master/src/python/motu-client.py", #path to 'motu-client.py' (https://github.com/clstoulouse/motu-client-python),
                            # Login Credentials
                            log_cmems=log,   
                            pwd_cmems=pass, 
                            # Motu Server and chosen Product/Dataset
                            motu_sc="http://nrtcmems.mercator-ocean.fr/mis-gateway-servlet/Motu",
                            serv_id="http://purl.org/myocean/ontology/service/database#GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS",
                            dataset_id="global-analysis-forecast-phy-001-024",
                            # Date 
                            yyyystart="2013",
                            mmstart="12",
                            # Area 
                            xmin=as.character(cellExt[1,1]),
                            xmax=as.character(cellExt[1,1]),
                            ymin=as.character(cellExt[1,2]),
                            ymax=as.character(cellExt[1,2]),
                            zsmall="0.494", 
                            zbig="1",
                            # Variables 
                            var = "thetao",
                            # Output files 
                            out_path = paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID,"/"), #Make sure to end your path with "/" 
                            pre_name= "monthly_"){
  
  
  #cell Extent
  cellExt <- getCellCentr(cellID)
  
  ## sourcing the function
  source("./Scripts/getCMEMS/getCMEMS.R")
  
  ## parameters
  
  ### motu-client.py path
  ### this may chnage according to your computer configuration
  motu_cl_lib <- "./Scripts/getCMEMS/libs/motu-client-python-master/src/python/motu-client.py"
  
  ### output dir
  # outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID,"/")
  
  ### credentials (in cred.txt)
  source("./Scripts/getCMEMS/cred.txt")
  
  
  
  
  getCMEMS_monthly(motu_cl = motu_cl_lib , 
                   out_path = out_path,
                   log_cmems = log_cmems,
                   pwd_cmems = pwd_cmems,
                    # Date 
                    yyyystart="2013",
                    mmstart="01",
                    # Area 
                    xmin=xmin,
                    xmax=xmax,
                    ymin=ymin,
                    ymax=ymax,
                    zsmall=zsmall, 
                    zbig=zbig,
                    #OutPath
                    out_path = out_path)
  
  
  
} # eo getCellTempData





#4. lapply the function that download temp data for a cell

TempData <- lapply(c(4622441,4657003), getCellTempData)

TempData <- lapply(4657003, getCellTempData)
 




res_monthly <- lapply(list(4622441, 4657003), getCMEMS_monthly(motu_cl = motu_cl_lib , 
                                out_path = outDir,
                                log_cmems = log,
                                pwd_cmems = pass,
                                # Date 
                                yyyystart="2013",
                                mmstart="01",
                                #extents
                                xmin="-176.6244",
                                xmax="-176.6244",
                                ymin="0.80875",
                                ymax="0.80875",
                                zsmall="0.494", 
                                zbig="1"))





TempData <- lapply(yearSeq,lapply(monthSeq, lapply(clist, function(cellID){
  
  #cell Extent
  cellExt <- getCellCentr(cellID)
  
  ## sourcing the function
  source("./Scripts/getCMEMS/getCMEMS.R")
  
  ## parameters
  
  ### motu-client.py path
  ### this may chnage according to your computer configuration
  motu_cl_lib <- "./Scripts/getCMEMS/libs/motu-client-python-master/src/python/motu-client.py"
  
  ### output dir
  # dir.create(path = paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID))
  # outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID,"/")
  outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/")
  
  ### credentials (in cred.txt)
  source("./Scripts/getCMEMS/cred.txt")
  
  prename= paste0("monthly_",cellID)
  
  
  res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib ,
                                  log_cmems = log,
                                  pwd_cmems = pass,
                                  # Date 
                                  yyyystart=yearSeq,
                                  mmstart=monthSeq,
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
  
  res_monthly
  
} # eo fct

)
)
)




















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
  dir.create(path = paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID))
  outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID,"/")
  
  ### credentials (in cred.txt)
  source("./Scripts/getCMEMS/cred.txt")
  
  
  
  
  res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib ,
                                  log_cmems = log,
                                  pwd_cmems = pass,
                                  # Date 
                                  yyyystart="2013",
                                  mmstart="01",
                                  # Area 
                                  xmin=as.character(cellExt[1,1]),
                                  xmax=as.character(cellExt[1,1]),
                                  ymin=as.character(cellExt[1,2]),
                                  ymax=as.character(cellExt[1,2]),
                                  zsmall="0.494", 
                                  zbig="1",
                                  #OutPath
                                  out_path = outDir)
  
  res_monthly
  
} # eo getCellTempData


test <- getCellTempData(4622441)


#4. lapply the function that download temp data for a cell

TempData <- lapply(list(4622441,4657003), getCellTempData)




#4. for loop for all cells

## sourcing the function
source("./Scripts/getCMEMS/getCMEMS.R")

### credentials (in cred.txt)
source("./Scripts/getCMEMS/cred.txt")

### motu-client.py path
### this may chnage according to your computer configuration
motu_cl_lib <- "./Scripts/getCMEMS/libs/motu-client-python-master/src/python/motu-client.py"

getCellTempData <- function(cellList){
  
  for(cellID in cellList){
    
    #cell Extent
    cellExt <- getCellCentr(cellID)
    
    # Output directory
    outDir <- paste0("/home/florian/NicheModelling_FBaletaud/data/rawdata/Environment/temp/CMEMS/",cellID,"/")
    
    #getCMEMS
    res_monthly[cellID] <- getCMEMS_monthly(motu_cl = motu_cl_lib ,
                                            log_cmems = log,
                                            pwd_cmems = pass,
                                            # Date 
                                            yyyystart="2013",
                                            mmstart="01",
                                            # Area 
                                            xmin=as.character(cellExt[1,1]),
                                            xmax=as.character(cellExt[1,1]),
                                            ymin=as.character(cellExt[1,2]),
                                            ymax=as.character(cellExt[1,2]),
                                            zsmall="0.494", 
                                            zbig="1",
                                            #OutPath
                                            out_path = outDir)
    
    
    
    
  }
  
}

tempData <- getCellTempData(cellList)



# reading cell files and storing them

cellFileList <- dir("./data/rawdata/Environment/temp/CMEMS")


#Function generating variable for one cell

# WARNING : This function might reach the limited number of files allowed to be opened at the same time 
# as the ncdf4 package doesnt seem to close the opened files during a loop function

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
    
    #Max
    maxTemp <- max(unlist(celldata))
    
    #Annual Range
    tempRange <- maxTemp - meanTemp
    
    
    df <- data.frame(cellID = cellID, MAT = meanTemp, MINT = minTemp, MAXT = maxTemp, TRAN = tempRange)
    
    # rownames(df) <- cellID
    
    # nc_close(cellFiles)
    
    # closeAllNcfiles()
    
    df
    
  }
  )
  
  
  
  data <- do.call("rbind", cellData)
  
  data
  
}


# tempVar <- getCellTempVar(clist = clist1, cellFileList)





############# Removing cells unneeded ###############



cellFileListToRemove <- as.vector(sapply(as.character(cellListAbs1),list.files,path = "./data/rawdata/Environment/temp/CMEMS", full.names = TRUE))

file.remove(cellFileListToRemove)

cellFileList <- dir("./data/rawdata/Environment/temp/CMEMS")


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









