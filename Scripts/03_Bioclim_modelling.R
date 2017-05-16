#################################################################
##                                                             ##
##         Modelling the first bioclimatic enveloppe           ##         
##                                                             ##
#################################################################

library(h2o)
library(h2oEnsemble)


earthGrid <- raster("./data/interdata/earthGrid.tif")
eezNcGrid <- raster("./data/interdata/eezNcGrid.tif")


## Reading datasets and binding with temperature cell data

spDatasets <- dir("./data/calibdata/bioclimodel", "*_speciesDataset.csv", full.names = TRUE)

spDatasets <- lapply(spDatasets, read.csv)

names(spDatasets) <- c("Carcharhinus_amblyrhynchos","Carcharhinus_melanopterus","Triaenodon_obesus")

cellTemp <- read.csv("./data/calibdata/bioclimodel/cellTemperatures.csv")
cellTemp <- cellTemp[,-1] # to be removed with the adding row.names = FALSE added in write.csv from occurrence_data.R
colnames(cellTemp)[1] <- "cellNumber" # To be removed after fixing scripts


dat <- lapply(spDatasets, function(tabs, preds = cellTemp) {
  
  tabs <- tabs[!duplicated(tabs$cellNumber), ]
  
  tabs$cellNumber <- cellFromXY(as.matrix(tabs[,c("decimalLongitude","decimalLatitude")]), object = earthGrid) # To remove with the fix between earthgrid cells and merged earthGrid/antiabs cells from occurrence_data.R
  
  preds <- preds[match(tabs$cellNumber, preds$cellNumber, nomatch=0),]
  
  df <- cbind(tabs, preds)
  
  df <- df[,c("occurrence","Tmin","Tmax","Tmean","Trange")]
  
  df$occurrence[df$occurrence == 1] <- "Presence"
  
  df$occurrence[df$occurrence == 0] <- "Absence"
  
  df <- df[sample(nrow(df), nrow(df)), ]
  
  df$occurrence <- factor(df$occurrence)
  
  df
  
  write.csv(df, file = paste0("./data/calibdata/bioclimodel/",tabs$species[1],"_calib.csv"), row.names = FALSE)
  
  cat(paste0("dataset written in"," ","./data/calibdata/bioclimodel/",tabs$species[1],"_calib.csv","\n"))
  
  return(df)
  
})

library(h2o)
library(h2oEnsemble)

# Initiating h2o

h2o.clust <- tryCatch(h2o.init(startH2O = FALSE),error=function(e){
  
  h2o.init(ip = "localhost", port = 54321, startH2O = TRUE,
           forceDL = FALSE, enable_assertions = TRUE, license = NULL,
           nthreads = 30, max_mem_size = NULL, min_mem_size = NULL,
           ice_root = tempdir(), strict_version_check = TRUE,
           proxy = NA_character_, https = FALSE, insecure = FALSE,
           cluster_id = NA_integer_, cookies = NA_character_)
  
  #h2o.removeAll() # (Optional) Remove all objects in H2O cluster
  
})


# Carcharhinus amblyrhynchos
datAmbly <- h2o.importFile("data/calibdata/bioclimodel/Carcharhinus_amblyrhynchos_calib.csv", destination_frame = "CamblyData", parse = TRUE, header = TRUE,
                      sep = ",", col.names = NULL, col.types = c("factor","numeric","numeric","numeric","numeric"), na.strings = NULL)



# Creating the 3 splits

splits <- h2o.splitFrame(
  datAmbly,           ##  splitting the H2O frame we read above
  c(0.7,0.15),   ##  create splits of 70% and 15%; 
  ##  H2O will create one more split of 1-(sum of these parameters)
  ##  so we will get 0.7 / 0.15 / 1 - (0.7+0.15) = 0.7/0.15/0.15
  seed=1234)    ##  setting a seed will ensure reproducible results (not R's seed)

train <- h2o.assign(splits[[1]], "train.hex")   
## assign the first result the R variable train
## and the H2O name train.hex
valid <- h2o.assign(splits[[2]], "valid.hex")   ## R valid, H2O valid.hex
test <- h2o.assign(splits[[3]], "test.hex")     ## R test, H2O test.hex

## take a look at the first few rows of the data set
train[1:5,]   ## rows 1-5, all columns




######### Ensemble Forecasting ###########

## Choosing the learners and metalearner
learner <- c("h2o.glm.wrapper", "h2o.randomForest.wrapper", 
             "h2o.gbm.wrapper", "h2o.deeplearning.wrapper") #
metalearner <- "h2o.glm.wrapper"


FitAmblyh2oEnsemble <- h2o.ensemble(x = names(train[-1]), y = "occurrence", 
                                        training_frame = train,
                                        family = "binomial",
                                        learner = learner,
                                        metalearner = metalearner,
                                        cvControl = list(V = 5))





## Saving results
h2o.save_ensemble(FitAmblyh2oEnsemble , path = "data/results/h2o_models/bioclimodel/Camblyrhynchos" , export_levelone = TRUE , force = TRUE)

### Evaluate Model Performance
perfAmbly <- h2o.ensemble_performance(FitAmblyh2oEnsemble, newdata = test)
perfAmbly
print(perfAmbly, metric = "MSE")







# Carcharhinus melanopterus

datMelano <- h2o.importFile("data/calibdata/bioclimodel/Carcharhinus_melanopterus_calib.csv", destination_frame = "CmelanoData", parse = TRUE, header = TRUE,
                      sep = ",", col.names = NULL, col.types = c("factor","numeric","numeric","numeric","numeric"), na.strings = NULL)

# Creating the 3 splits

splits <- h2o.splitFrame(
  datMelano,           ##  splitting the H2O frame we read above
  c(0.7,0.15),   ##  create splits of 70% and 15%; 
  ##  H2O will create one more split of 1-(sum of these parameters)
  ##  so we will get 0.7 / 0.15 / 1 - (0.7+0.15) = 0.7/0.15/0.15
  seed=1234)    ##  setting a seed will ensure reproducible results (not R's seed)

train <- h2o.assign(splits[[1]], "train.hex")   
## assign the first result the R variable train
## and the H2O name train.hex
valid <- h2o.assign(splits[[2]], "valid.hex")   ## R valid, H2O valid.hex
test <- h2o.assign(splits[[3]], "test.hex")     ## R test, H2O test.hex

## take a look at the first few rows of the data set
train[1:5,]   ## rows 1-5, all columns




######### Ensemble Forecasting ###########


FitMelanoh2oEnsemble <- h2o.ensemble(x = names(train[-1]), y = "occurrence", 
                    training_frame = train,
                    family = "binomial",
                    learner = learner,
                    metalearner = metalearner,
                    cvControl = list(V = 5))





## Saving results
h2o.save_ensemble(FitMelanoh2oEnsemble , path = "data/results/h2o_models/bioclimodel/Cmelanopterus" , export_levelone = TRUE , force = TRUE)

### Evaluate Model Performance
perfMelano <- h2o.ensemble_performance(FitMelanoh2oEnsemble, newdata = test)
perfMelano
print(perfMelano, metric = "MSE")




# Triaenodon obesus

datObesus <- h2o.importFile("data/calibdata/bioclimodel/Triaenodon_obesus_calib.csv", destination_frame = "TobesusData", parse = TRUE, header = TRUE,
                            sep = ",", col.names = NULL, col.types = c("factor","numeric","numeric","numeric","numeric"), na.strings = NULL)

# Creating the 3 splits

splits <- h2o.splitFrame(
  datObesus,           ##  splitting the H2O frame we read above
  c(0.7,0.15),   ##  create splits of 70% and 15%; 
  ##  H2O will create one more split of 1-(sum of these parameters)
  ##  so we will get 0.7 / 0.15 / 1 - (0.7+0.15) = 0.7/0.15/0.15
  seed=1234)    ##  setting a seed will ensure reproducible results (not R's seed)

train <- h2o.assign(splits[[1]], "train.hex")   
## assign the first result the R variable train
## and the H2O name train.hex
valid <- h2o.assign(splits[[2]], "valid.hex")   ## R valid, H2O valid.hex
test <- h2o.assign(splits[[3]], "test.hex")     ## R test, H2O test.hex

## take a look at the first few rows of the data set
train[1:5,]   ## rows 1-5, all columns




######### Ensemble Forecasting ###########


FitObesush2oEnsemble <- h2o.ensemble(x = names(train[-1]), y = "occurrence", 
                         training_frame = train,
                         family = "binomial",
                         learner = learner,
                         metalearner = metalearner,
                         cvControl = list(V = 5))





## Saving results
h2o.save_ensemble(FitObesush2oEnsemble , path = "data/results/h2o_models/bioclimodel/TObesus" , export_levelone = TRUE , force = TRUE)

### Evaluate Model Performance
perfObesus <- h2o.ensemble_performance(FitObesush2oEnsemble, newdata = test)
perfObesus
print(perfObesus, metric = "MSE")




############ Predicting DATA #############

FitAmblyh2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/bioclimodel/Camblyrhynchos")



# Amblyrhynchos
predictData <- stack(list.files("./data/predictdata/BEM", pattern = ".tif", full.names = TRUE))
predData.R <- as.data.frame(predictData)
predData.R <- na.omit(predData.R)
write.csv(predData.R, "./data/predictdata/BEM/predDataBEM.csv", row.names = FALSE)
predDatah2o <- h2o.importFile("./data/predictdata//BEM/predDataBEM.csv", destination_frame = "predDataBEM", parse = TRUE, header = TRUE,
                              sep = ",", col.names = NULL, col.types = c("numeric","numeric","numeric","numeric"), na.strings = NULL)

amblyPreds <- predict(FitAmblyh2oEnsemble, predDatah2o)


predData.R <- as.data.frame(predictData)
predData.R <- cbind(as.data.frame(coordinates(predictData)), predData.R)
predDataCoords <- na.omit(predData.R)
predDataCoords <- as.matrix(predDataCoords[,c(1:2)])

predictions <- as.data.frame(amblyPreds$pred)

amblyMapPoints <- SpatialPointsDataFrame(
  coords = predDataCoords,
  data = predictions,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



amblyMapGrid <- amblyMapPoints
gridded(amblyMapGrid) = TRUE
amblyMapGrid <- as(amblyMapGrid, "SpatialGridDataFrame")
amblyMap <- raster(amblyMapGrid, values = TRUE)
amblyMap <- rasterize(amblyMapPoints, amblyMap, field = "Presence")


plot(amblyMap)
writeRaster(amblyMap, filename = paste0("./data/predictdata/Filtered/BEM_ambly.tif"), overwrite = TRUE)




# Melanopterus

FitMelanoh2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/bioclimodel/Cmelanopterus")

predictData <- stack(list.files("./data/predictdata/BEM", pattern = ".tif", full.names = TRUE))
predData.R <- as.data.frame(predictData)
predData.R <- na.omit(predData.R)
write.csv(predData.R, "./data/predictdata/BEM/predDataBEM.csv", row.names = FALSE)
predDatah2o <- h2o.importFile("./data/predictdata//BEM/predDataBEM.csv", destination_frame = "predDataBEM", parse = TRUE, header = TRUE,
                              sep = ",", col.names = NULL, col.types = c("numeric","numeric","numeric","numeric"), na.strings = NULL)

melanoPreds <- predict(FitMelanoh2oEnsemble, predDatah2o)


predData.R <- as.data.frame(predictData)
predData.R <- cbind(as.data.frame(coordinates(predictData)), predData.R)
predDataCoords <- na.omit(predData.R)
predDataCoords <- as.matrix(predDataCoords[,c(1:2)])

predictions <- as.data.frame(melanoPreds$pred)

melanoMapPoints <- SpatialPointsDataFrame(
  coords = predDataCoords,
  data = predictions,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



melanoMapGrid <- melanoMapPoints
gridded(melanoMapGrid) = TRUE
melanoMapGrid <- as(melanoMapGrid, "SpatialGridDataFrame")
melanoMap <- raster(melanoMapGrid, values = TRUE)
melanoMap <- rasterize(melanoMapPoints, melanoMap, field = "Presence")


plot(melanoMap)
writeRaster(melanoMap, filename = paste0("./data/predictdata/Filtered/BEM_melano.tif"), overwrite = TRUE)


# Obesus

FitObesush2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/bioclimodel/TObesus")

predictData <- stack(list.files("./data/predictdata/BEM", pattern = ".tif", full.names = TRUE))
predData.R <- as.data.frame(predictData)
predData.R <- na.omit(predData.R)
write.csv(predData.R, "./data/predictdata/BEM/predDataBEM.csv", row.names = FALSE)
predDatah2o <- h2o.importFile("./data/predictdata//BEM/predDataBEM.csv", destination_frame = "predDataBEM", parse = TRUE, header = TRUE,
                              sep = ",", col.names = NULL, col.types = c("numeric","numeric","numeric","numeric"), na.strings = NULL)

obesusPreds <- predict(FitObesush2oEnsemble, predDatah2o)


predData.R <- as.data.frame(predictData)
predData.R <- cbind(as.data.frame(coordinates(predictData)), predData.R)
predDataCoords <- na.omit(predData.R)
predDataCoords <- as.matrix(predDataCoords[,c(1:2)])

predictions <- as.data.frame(obesusPreds$pred)

obesusMapPoints <- SpatialPointsDataFrame(
  coords = predDataCoords,
  data = predictions,
  proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))



obesusMapGrid <- obesusMapPoints
gridded(obesusMapGrid) = TRUE
obesusMapGrid <- as(obesusMapGrid, "SpatialGridDataFrame")
obesusMap <- raster(obesusMapGrid, values = TRUE)
obesusMap <- rasterize(obesusMapPoints, obesusMap, field = "Presence")


plot(obesusMap)
writeRaster(obesusMap, filename = paste0("./data/predictdata/Filtered/BEM_obesus.tif"), overwrite = TRUE)


# Shut down H2o Cluster
h2o.shutdown()
Y
