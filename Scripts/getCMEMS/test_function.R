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
outDir <- "./Scripts/getCMEMS/downs/"

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
	    mmend="03",
	    hh=" 12:00:00",
	    dd="01",
	    xmin="163.473146",
	    xmax="163.633146",
	    ymin="-21.523697",
	    ymax="-21.363697")



res_monthly <- getCMEMS_monthly(motu_cl = motu_cl_lib , 
		        out_path = outDir,
		        log_cmems = log,
		        pwd_cmems = pass,
		        # Date 
		        yyyystart="2013",
		        mmstart="01")

## ---- ReadResults
library(ncdf4)


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

#Read the time
time_mon <- ncvar_get(mon, "time")

length(time_mon) == length(res_monthly)

