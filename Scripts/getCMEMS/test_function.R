################################################
# Testing :
# R wrapper for motu CMEMS data download
#
# francois.guilhaumon@ird.fr
################################################

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

### call the functions
res <- getCMEMS(motu_cl = motu_cl_lib , 
	    out_path = outDir,
	    log_cmems = log,
	    pwd_cmems = pass,
	    # Date 
	    yyyystart="2013",
	    mmstart="01",
	    yyyyend="2013",
	    mmend="02",
	    hh=" 12:00:00",
	    dd="01",
	    xmin="-176.6244",
	    xmax="-176.6244",
	    ymin="0.80875",
	    ymax="0.80875",
	    zsmall="0.494", 
	    zbig="1")



res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib , 
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
		        zbig="1")

## ---- ReadResults
library(ncdf4)


#getCMEMS
fev <- nc_open("./data/rawdata/Environment/temp/CMEMS/monthly_3048157global-analysis-forecast-phy-001-024_thetao_2013-02.nc")
fev2 <- nc_open("./data/rawdata/Environment/temp/CMEMS/monthly_3201152global-analysis-forecast-phy-001-024_thetao_2016-02.nc")

# Read longitude & latitude
lon <- ncvar_get(jan, "longitude")
lat <- ncvar_get(jan, "latitude")

#Read the time
time_jan <- ncvar_get(jan, "time")

length(time_jan) == length(res)

# Read values
temp_fev <- ncvar_get(fev,"thetao")
temp_fev2 <- ncvar_get(fev2,"thetao")

#monthly

mon <- nc_open("./Scripts/getCMEMS/downs/monthly_global-analysis-forecast-phy-001-024_thetao_2013-01.nc")

#Read the time
time_mon <- ncvar_get(mon, "time")

length(time_mon) == length(res_monthly)

