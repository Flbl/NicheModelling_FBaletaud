#################################################################
##                                                             ##
##         Modelling the second regional enveloppe             ##         
##                                                             ##
#################################################################

library(h2o)
library(h2oEnsemble)


## Reading datasets and binding with temperature cell data

spDatasets <- dir("./data/calibdata/regionmodel", "*_speciesDataset_region.csv", full.names = TRUE)

spDatasets <- lapply(spDatasets, read.csv)

names(spDatasets) <- c("Carcharhinus_amblyrhynchos","Carcharhinus_melanopterus","Triaenodon_obesus")

dat <- lapply(spDatasets, function(tabs) {
  
  df <- tabs

  df$occurrence[df$occurrence == 1] <- "Presence"

  df$occurrence[df$occurrence == 0] <- "Absence"

  df$occurrence <- factor(df$occurrence)
  
  drops <- c("species","decimalLongitude","decimalLatitude","cellNumber")
  
  df <- df[, !(names(df) %in% drops)]
  
  df <- df[sample(nrow(df), nrow(df)), ]

  write.csv(df, file = paste0("./data/calibdata/regionmodel/",tabs$species[1],"_region_calib.csv"), row.names = FALSE)
  
  cat(paste0("dataset written in"," ","./data/calibdata/regionmodel/",df$species[1],"_region_calib.csv","\n"))
  
  return(df)
  
})


## Launching h2o


h2o.init()


# Carcharhinus amblyrhynchos
datAmbly <- h2o.importFile("data/calibdata/regionmodel/Carcharhinus_amblyrhynchos_region_calib.csv", destination_frame = "CamblyData", parse = TRUE, header = TRUE,
                           sep = ",", col.names = NULL, col.types = c("factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"), na.strings = NULL)



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


ESAmbly <- h2o.ensemble(x = names(train[-1]), y = "occurrence", 
                        training_frame = train,
                        family = "binomial",
                        learner = learner,
                        metalearner = metalearner,
                        cvControl = list(V = 5))





## Saving results
h2o.save_ensemble(ESAmbly , path = "data/results/h2o_models/Camblyrhynchos/" , export_levelone = TRUE , force = TRUE)

### Evaluate Model Performance
perfAmbly <- h2o.ensemble_performance(ESAmbly, newdata = test)
perfAmbly
print(perfAmbly, metric = "MSE")











# Carcharhinus melanopterus # DATA SET TOO SMALL

datMelano <- h2o.importFile("data/calibdata/regionmodel/Carcharhinus_melanopterus_region_calib.csv", destination_frame = "CmelanoData", parse = TRUE, header = TRUE,
                            sep = ",", col.names = NULL, col.types = c("factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"), na.strings = NULL)

# Creating the 3 splits

splits <- h2o.splitFrame(
  datMelano,           ##  splitting the H2O frame we read above
  c(0.9,0.05),   ##  create splits of 70% and 15%;
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


ESMelano <- h2o.ensemble(x = names(train[-1]), y = "occurrence",
                         training_frame = train,
                         family = "binomial",
                         learner = learner,
                         metalearner = metalearner,
                         cvControl = list(V = 5))





## Saving results
h2o.save_ensemble(ESMelano , path = "data/results/h2o_models/Cmelanopterus" , export_levelone = TRUE , force = TRUE)

### Evaluate Model Performance
perfMelano <- h2o.ensemble_performance(ESMelano, newdata = test)
perfMelano
print(perfMelano, metric = "MSE")


# Triaenodon obesus

datObesus <- h2o.importFile("data/calibdata/regionmodel/Triaenodon_obesus_region_calib.csv", destination_frame = "TobesusData", parse = TRUE, header = TRUE,
                            sep = ",", col.names = NULL, col.types = c("factor","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"), na.strings = NULL)

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


ESObesus <- h2o.ensemble(x = names(train[-1]), y = "occurrence", 
                         training_frame = train,
                         family = "binomial",
                         learner = learner,
                         metalearner = metalearner,
                         cvControl = list(V = 5))





## Saving results
h2o.save_ensemble(ESObesus , path = "data/results/h2o_models/TObesus" , export_levelone = TRUE , force = TRUE)

### Evaluate Model Performance
perfObesus <- h2o.ensemble_performance(ESObesus, newdata = test)
perfObesus
print(perfObesus, metric = "MSE")


##### PREDICTING DATA #####


predictData <- stack(list.files("./data/predictdata", full.names = TRUE))
predData.R <- as.data.frame(predictData)
predData.R <- na.omit(predData.R)
write.csv(predData.R, "./data/predictdata/predDataRegion.csv", row.names = FALSE)
predDatah2o <- h2o.importFile("./data/predictdata/predDataRegion.csv", destination_frame = "predDataRegion", parse = TRUE, header = TRUE,
                                sep = ",", col.names = NULL, col.types = c("numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"), na.strings = NULL)

amblyPreds <- predict(ESAmbly, predDatah2o)



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

