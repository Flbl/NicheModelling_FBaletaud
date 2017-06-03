# Get results

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



FitBEMAmblyh2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/bioclimodel/Camblyrhynchos")

FitHABamblyh2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/regionmodel/Camblyrhynchos")


FitBEMmelanoh2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/bioclimodel/Cmelanopterus")

FitHABmelanoh2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/regionmodel/Cmelanopterus")

FitBEMobesush2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/bioclimodel/TObesus")

FitHABobesush2oEnsemble <- h2o.load_ensemble("./data/results/h2o_models/regionmodel/TObesus")



mod = FitHABamblyh2oEnsemble

getResults <- function(mod){
  
  varImpList <- sapply(mod$basefits, h2o.varimp)
  
  varImpList$h2o.glm.wrapper <- varImpList$h2o.glm.wrapper[,c(1,2)]
  varImpList$h2o.glm.wrapper <- na.omit(varImpList$h2o.glm.wrapper)
  varImpList$h2o.randomForest.wrapper <- varImpList$h2o.randomForest.wrapper[, c(1,3)]
  varImpList$h2o.gbm.wrapper <- varImpList$h2o.gbm.wrapper[,c(1,3)]
  varImpList$h2o.deeplearning.wrapper <- NULL
  
  
  par(mfrow = c(2,2))
  par(mar=c(4,10,2,2))
  
  
  barplot(rev(varImpList$h2o.glm.wrapper[,2]), names.arg = rev(varImpList$h2o.glm.wrapper[,1]), horiz = TRUE, col = "dodgerblue4", las = 1, xlab = "Importance des variables (coefficients)" )
  title(main = "A", adj = 0)
  barplot(rev(varImpList$h2o.randomForest.wrapper[,2]), names.arg = rev(varImpList$h2o.randomForest.wrapper[,1]), horiz = TRUE, col = "dodgerblue4", las = 1, xlab = "Importance des variables (réduites)")
  title(main = "B", adj = 0)
  barplot(rev(varImpList$h2o.gbm.wrapper[,2]), names.arg = rev(varImpList$h2o.gbm.wrapper[,1]), horiz = TRUE, col = "dodgerblue4", las = 1, xlab = "Importance des variables (réduites)")
  title( main = "C", adj = 0)
  
  
  
  # lapply(varImpList, function(x) {
  # 
  #   barplot(rev(x[,2]), names.arg = rev(x[,1]), horiz = TRUE, col = "dodgerblue4", las = 1, xlab = "Importance des variables")
  # 
  # })
  
  
  
  modImp <- h2o.varimp(mod$metafit)
  modImp <- na.omit(modImp)
  
  barplot(rev(modImp[, 2]), names.arg = rev(modImp[,1]), horiz = TRUE, col = "dodgerblue4", las = 1, xlab = "Importance des modèles (coefficients)")
  title( main = "D", adj = 0)
  
  
  
  
  confmatList <- lapply(mod$basefits, h2o.confusionMatrix)
  metrcs <- lapply(confmatList, function(cfm) {
    
    mtrcs <- 1 - cfm$Error
    
    names(mtrcs) <- c("Specificité","Sensitivité","Précision")
    TSS <- (mtrcs["Sensitivité"] + mtrcs["Specificité"]) - 1
    names(TSS) <- "TSS"
    mtrcs <- append(mtrcs, TSS)
    mtrcs
    
  })
  
  metrcsES <- h2o.confusionMatrix(mod$metafit)
  metrcsES <- 1 - metrcsES$Error
  names(metrcsES) <- c("Specificité","Sensitivité","Précision")
  TSS <- (metrcsES["Sensitivité"] + metrcsES["Specificité"]) - 1
  names(TSS) <- "TSS"
  
  metrcsES <- append(metrcsES, TSS)
  metrcsES <- list(metrcsES)
  names(metrcsES) <- "Ensemble Stacking"
  
  c(metrcs, metrcsES)
  
  

}

getResults(FitBEMAmblyh2oEnsemble)

getResults(HABamblyFit)

getResults(FitBEMmelanoh2oEnsemble)

getResults(FitHABmelanoh2oEnsemble)

getResults(FitBEMobesush2oEnsemble)

getResults(FitHABobesush2oEnsemble)


getGlobVarImp <- function(mod){
  
  varImpList <- sapply(mod$basefits, h2o.varimp)
  varImpList$h2o.deeplearning.wrapper <- NULL
  varImpList <- lapply(varImpList, as.data.frame)
  varImpList <- lapply(varImpList,na.omit)
  
  varImpList$h2o.glm.wrapper$percentage <- sapply(varImpList$h2o.glm.wrapper$coefficients,function(x) x/sum(varImpList$h2o.glm.wrapper$coefficients))
  
  # same variable order
  varImpList <- lapply(varImpList, function(tab) {
    
    rownames(tab) <- tab[,1]
    
    if("Tmax" %in% rownames(tab))  tab <- tab[c("Tmax","Tmean","Tmin","Trange"),]
    
    if("Abyss" %in% rownames(tab)) tab <- tab[c("Coral","Abyss","Travel_Dist","Slope","Plateaus","Canyons","Shelf","Ridges","Basins","Escarpments","Guyots"),]

    tab
    
  })
  
  # dataframe with the different models
  
  varImpdf <- data.frame(variable = varImpList[[1]]$names, 
                         h2o.glm.wrapper = varImpList[[1]]$percentage, 
                         h2o.randomForest.wrapper = varImpList[[2]]$percentage, 
                         h2o.gbm.wrapper = varImpList[[3]]$percentage)
  
  modImp <- h2o.varimp(mod$metafit)
  modImp <- na.omit(modImp)
  modImp <- as.data.frame(modImp)
  modImp <- modImp[-which(modImp$names == "h2o.deeplearning.wrapper"),]
  modImp$modelCoeffs <- sapply(modImp$coefficients,function(x) x/sum(modImp$coefficients))
  
  modImpdf <- data.frame(model = modImp$names, coeffs = modImp$modelCoeffs)
  
  varImpdf$weightedGLM <- sapply(varImpdf$h2o.glm.wrapper, function(x) x = x*modImpdf[which(modImpdf$model == "h2o.glm.wrapper"), "coeffs"])
  varImpdf$weightedDRF <- sapply(varImpdf$h2o.randomForest.wrapper, function(x) x = x*modImpdf[which(modImpdf$model == "h2o.randomForest.wrapper"), "coeffs"])
  varImpdf$weightedGBM <- sapply(varImpdf$h2o.gbm.wrapper, function(x) x = x*modImpdf[which(modImpdf$model == "h2o.gbm.wrapper"), "coeffs"])
  
  varImpdf$weightedVarImp <- rowSums(varImpdf[,c("weightedGLM","weightedDRF","weightedGBM")])
  
  barplot(varImpdf$weightedVarImp , names.arg = varImpdf$variable,cex.names = 1.3, ylab = "importance des variables (réduites)", col = "orangered3", las = 2, xpd= FALSE, ylim = c(0,0.5))
  # text(seq(length(varImpdf$weightedVarImp)),par("usr")[3], labels = varImpdf$variable, srt = 45, adj = 1, cex = 1.2,
  #      xpd = TRUE)
}

par(mar = c(9, 4, 1, 2))
par(mfcol = c(2,3))

getGlobVarImp(FitBEMAmblyh2oEnsemble)

getGlobVarImp(FitHABamblyh2oEnsemble)

getGlobVarImp(FitBEMmelanoh2oEnsemble)

getGlobVarImp(FitHABmelanoh2oEnsemble)

getGlobVarImp(FitBEMobesush2oEnsemble)

getGlobVarImp(FitHABobesush2oEnsemble)






# Shut down H2o Cluster
h2o.shutdown()
Y

