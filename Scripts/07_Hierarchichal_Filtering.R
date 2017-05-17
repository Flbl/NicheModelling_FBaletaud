#################################################################
##                                                             ##
##      Multiplying enveloppes to get final predictions        ##         
##                                                             ##
#################################################################



eezNcPoly <- readOGR(dsn = "./data/interdata", layer = "eezNcPoly")

tifs <- list.files("./data/predictdata/Filtered", full.names = TRUE)

tifs <- stack(list.files("./data/predictdata/Filtered", full.names = TRUE))

# Amblyrhynchos

plot(tifs[[1]])
plot(eezNcPoly, add = TRUE)

plot(tifs[[4]])
plot(eezNcPoly, add = TRUE)

filtAmblyMap <- tifs[[1]] * tifs[[4]]
plot(filtAmblyMap)
plot(eezNcPoly, add = TRUE)


# Melanopterus

plot(tifs[[2]])
plot(eezNcPoly, add = TRUE)

plot(tifs[[5]])
plot(eezNcPoly, add = TRUE)

filtmelanoMap <- tifs[[2]] * tifs[[5]]
plot(filtmelanoMap)
plot(eezNcPoly, add = TRUE)



# Melanopterus

plot(tifs[[3]])
plot(eezNcPoly, add = TRUE)

plot(tifs[[6]])
plot(eezNcPoly, add = TRUE)

filtObesusMap <- tifs[[3]] * tifs[[6]]
plot(filtObesusMap)
plot(eezNcPoly, add = TRUE)
