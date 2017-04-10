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

nc = nc_open(filename = "./Environment/temp/rawData/4622441/CMEMS_GLO_001_024_TEMP_2007-01.nc")
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



############### Long Method that needs to be changed ####################
clist1.4 <- clist[1:4]

clist5.8 <- clist[5:8]

clist9.12 <- clist[9:12]

clist13.16 <- clist[13:16]

clist17.20 <- clist[17:20]

clist21.24 <- clist[21:24]

clist25.28 <- clist[25:28]

clist29.32 <- clist[29:32]

clist33.36 <- clist[33:36]

clist37.40 <- clist[37:40]

clist41.44 <- clist[41:44]

clist45.48 <- clist[45:48]

clist49.52 <- clist[49:52]

clist53.56 <- clist[53:56]

clist57.60 <- clist[57:60]

clist61.64 <- clist[61:64]

clist65.68 <- clist[65:68]

clist69.72 <- clist[69:72]

clist73.76 <- clist[73:76]

clist77.80 <- clist[77:80]

clist81.84 <- clist[81:84]

clist85.88 <- clist[85:88]

clist89.92 <- clist[89:92]

clist93.96 <- clist[93:96]

clist97.100 <- clist[97:100]

clist101.104 <- clist[101:104]

clist105.108 <- clist[105:108]

clist109.112 <- clist[109:112]

clist113.116 <- clist[113:116]

clist117.120 <- clist[117:120]

clist121.124 <- clist[121:124]

clist125.128 <- clist[125:128]

clist129.132 <- clist[129:132]

clist133.136 <- clist[133:136]

clist137.140 <- clist[137:140]

clist141.142 <- clist[141:142]



tempVar4 <- getCellTempVar(clist = clist1.4, cellFileList)
tempVar8 <- getCellTempVar(clist = clist5.8, cellFileList)
tempVar12 <- getCellTempVar(clist = clist9.12, cellFileList)
tempVar16 <- getCellTempVar(clist = clist13.16, cellFileList)
tempVar20 <- getCellTempVar(clist = clist17.20, cellFileList)
tempVar24 <- getCellTempVar(clist = clist21.24, cellFileList)
tempVar28 <- getCellTempVar(clist = clist25.28, cellFileList)
tempVar32 <- getCellTempVar(clist = clist29.32, cellFileList)
tempVar36 <- getCellTempVar(clist = clist33.36, cellFileList)
tempVar40 <- getCellTempVar(clist = clist37.40, cellFileList)
tempVar44 <- getCellTempVar(clist = clist41.44, cellFileList)
tempVar48 <- getCellTempVar(clist = clist45.48, cellFileList)
tempVar52 <- getCellTempVar(clist = clist49.52, cellFileList)
tempVar56 <- getCellTempVar(clist = clist53.56, cellFileList)
tempVar60 <- getCellTempVar(clist = clist57.60, cellFileList)
tempVar64 <- getCellTempVar(clist = clist61.64, cellFileList)
tempVar68 <- getCellTempVar(clist = clist65.68, cellFileList)
tempVar72 <- getCellTempVar(clist = clist69.72, cellFileList)
tempVar76 <- getCellTempVar(clist = clist73.76, cellFileList)
tempVar80 <- getCellTempVar(clist = clist77.80, cellFileList)
tempVar84 <- getCellTempVar(clist = clist81.84, cellFileList)
tempVar88 <- getCellTempVar(clist = clist85.88, cellFileList)
tempVar92 <- getCellTempVar(clist = clist89.92, cellFileList)
tempVar96 <- getCellTempVar(clist = clist93.96, cellFileList)
tempVar100 <- getCellTempVar(clist = clist97.100, cellFileList)
tempVar104 <- getCellTempVar(clist = clist101.104, cellFileList)
tempVar108 <- getCellTempVar(clist = clist105.108, cellFileList)
tempVar112 <- getCellTempVar(clist = clist109.112, cellFileList)
tempVar116 <- getCellTempVar(clist = clist113.116, cellFileList)
tempVar120 <- getCellTempVar(clist = clist117.120, cellFileList)
tempVar124 <- getCellTempVar(clist = clist121.124, cellFileList)
tempVar128 <- getCellTempVar(clist = clist125.128, cellFileList)
tempVar132 <- getCellTempVar(clist = clist129.132, cellFileList)
tempVar136 <- getCellTempVar(clist = clist133.136, cellFileList)
tempVar140 <- getCellTempVar(clist = clist137.140, cellFileList)
tempVar142 <- getCellTempVar(clist = clist141.142, cellFileList)


tempVar <- rbind(tempVar4,tempVar8,tempVar12,tempVar16,tempVar20,tempVar24,tempVar28,
                 tempVar32,tempVar36,tempVar40,tempVar44,tempVar48,tempVar52,tempVar56,
                 tempVar60,tempVar64,tempVar68,tempVar72,tempVar76,tempVar80,tempVar84,
                 tempVar88,tempVar92,tempVar96,tempVar100,tempVar104,tempVar108,tempVar112,
                 tempVar116,tempVar120,tempVar124,tempVar128,tempVar132,tempVar136,tempVar140,tempVar142)


############### eo Long Method that needs to be changed ####################













