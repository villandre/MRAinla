## code to prepare `MODISdataTraining` and `MODISdataTest` dataset goes here

########## Function definitions ###############

uniformiseLandCover <- function(landCoverPointsList) {
  landCovers <- do.call("c", lapply(landCoverPointsList, function(x) colnames(x@data)))
  uniqueLandCovers <- unique(landCovers)
  landCoverIndices <- as.numeric(substr(uniqueLandCovers, start = 10, stop = 100))
  uniqueLandCovers <- uniqueLandCovers[order(landCoverIndices)]
  lapply(landCoverPointsList, function(landCoverPoints) {
    if (length(missingCols <- setdiff(uniqueLandCovers, colnames(landCoverPoints@data))) > 0) {
      landCoverPoints@data[missingCols] <- 0
      landCoverPoints@data <- landCoverPoints@data[ , uniqueLandCovers]
    }
    landCoverPoints
  })
}

prepareDataForMRAinla <- function(temperatures, elevations, landCover, satelliteNamesVec, collectionDatesPOSIX, completeDateVector = collectionDatesPOSIX) {
  if ("RasterLayer" %in% class(temperatures[[1]])) {
    temperaturePoints <- lapply(temperatures, FUN = raster::rasterToPoints, spatial = TRUE)
  } else {
    temperaturePoints <- temperatures
  }

  satelliteNamesList <- lapply(seq_along(satelliteNamesVec), function(dayIndex) {
    rep(satelliteNamesVec[[dayIndex]], nrow(temperaturePoints[[dayIndex]]@coords))
  })
  satellite <- do.call("c", satelliteNamesList)
  satellite <- as.numeric(factor(x = satellite, levels = c("Terra", "Aqua"))) - 1
  numTimePoints <- length(completeDateVector)

  timeValues <- do.call("c", lapply(seq_along(collectionDatesPOSIX), function(x) rep(collectionDatesPOSIX[[x]], length(temperaturePoints[[x]]))))

  timeLevels <- as.numeric(factor(as.character(timeValues), levels = as.character(completeDateVector))) - 1

  timeModelMatrix <- t(sapply(timeLevels, function(x) {
    unitVector <- rep(0, numTimePoints - 1)
    unitVector[x] <- 1 # When x takes value 0, the vector remains all 0s, which is what we want.
    unitVector
  }))
  colnames(timeModelMatrix) <- paste("time", 2:numTimePoints, sep = "")

  funToGetLandCoverPoints <- function(tempPoints) {
    tempPointsReproj <- sp::spTransform(tempPoints, crs(landCover))
    landCoverAtPoints <- raster::extract(landCover, tempPointsReproj)
    landCoverValues <- sort(unique(landCoverAtPoints))
    columnNames <- paste("landCover", landCoverValues, sep = "")
    landCoverMatrix <- t(sapply(landCoverAtPoints, function(x) {
      unitVec <- numeric(length(columnNames))
      unitVec[match(x, landCoverValues)] <- 1
      unitVec
    }))
    colnames(landCoverMatrix) <- columnNames
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = as.data.frame(landCoverMatrix), proj4string = raster::crs(tempPoints)) # A line in data with only zeros corresponds to a missing value.
  }

  landCoverPoints <- lapply(temperaturePoints, FUN = funToGetLandCoverPoints)
  landCoverPoints <- uniformiseLandCover(landCoverPoints)

  elevationPoints <- lapply(temperaturePoints, function(tempPoints) {
    tempPoints <- sp::spTransform(tempPoints, crs(elevations[[1]]))
    elevationValues <- rep(0, length(tempPoints))
    lapply(elevations, function(elevationRaster) {
      extractedValues <- raster::extract(elevationRaster, tempPoints)
      elevationValues[!is.na(extractedValues)] <<- extractedValues[!is.na(extractedValues)]
      NULL
    })
    sp::SpatialPointsDataFrame(coords = tempPoints@coords, data = data.frame(elevation = elevationValues), proj4string = crs(tempPoints))
  })

  combinedData <- do.call("cbind", lapply(list(landCoverPoints, elevationPoints), function(x) do.call("rbind", lapply(x, function(y) y@data))))
  combinedData <- cbind(combinedData, timeModelMatrix, Aqua = satellite)
  if ("RasterLayer" %in% class(temperatures[[1]])) {
    combinedData <- cbind(do.call("rbind", lapply(temperaturePoints, function(y) y@data)), combinedData)
    colnames(combinedData)[[1]] <- "y"
  }

  coordinates <- do.call("rbind", lapply(temperaturePoints, function(x) x@coords))
  rownames(coordinates) <- as.character(1:nrow(coordinates))
  missingLandCoverOrElevation <- (rowSums(combinedData[ , grep(colnames(combinedData), pattern = "landCover", value = TRUE)]) == 0) | is.na(combinedData[, "elevation"])

  spacetime::STIDF(sp = SpatialPoints(coordinates[!missingLandCoverOrElevation, ], proj4string = crs(temperaturePoints[[1]])), time = timeValues[!missingLandCoverOrElevation], data = as.data.frame(combinedData[!missingLandCoverOrElevation,]))
}

funToCreateRaster <- function(temperatureSdsList, polygonBound) {
  extractionFun <- function(x) {
    tempGrid <- readGDAL(x$SDS4gdal[1], as.is = TRUE)
    hourGrid <- readGDAL(x$SDS4gdal[3], as.is = TRUE)
    tempGrid$band1 <- tempGrid$band1 * 0.02 - 273.15 # See https://gis.stackexchange.com/questions/72524/how-do-i-convert-the-lst-values-on-the-modis-lst-image-to-degree-celsius
    # There's a 0.02 scaling factor applied to values in file to get the temperatures.
    # The -273.15 brings temperatures back in Celsius
    hourGrid@data[,1] <- hourGrid@data[,1] * 0.1

    list(temperatureRaster = raster(tempGrid), hourRaster = raster(hourGrid))
  }
  tempAndTimeRasters <- lapply(temperatureSdsList, extractionFun)

  createRaster <- function(rasterName) {
    rasterList <- lapply(tempAndTimeRasters, function(x) x[[rasterName]])
    mergedRasters <- do.call(raster::merge, rasterList)
    smallerRaster <- crop(x = mergedRasters, y = polygonBound)
    spObject <- rasterToPoints(smallerRaster, spatial = TRUE)
    polygonValuesIndex <- over(x = spObject, y = polygonBound)
    pointsInPolygon <- subset(spObject, subset = !is.na(polygonValuesIndex))
    values(smallerRaster) <- rep(NA, ncell(smallerRaster))
    if (sum(!is.na(polygonValuesIndex)) == 0) {
      return(smallerRaster)
    }
    rasterize(x = pointsInPolygon, y = smallerRaster, field = "layer")
  }
  rasterNames <- c("temperatureRaster", "hourRaster")
  tempAndTime <- lapply(rasterNames, FUN = createRaster)
  names(tempAndTime) <- rasterNames
  tempAndTime
}

#######################################################

library(sp)
library(raster)
library(MODIS)
library(rgdal)
library(rgeos)
library(spacetime)
library(RhpcBLASctl)
library(geoR)
library(maptools)

RandomFields::RFoptions(cores = 1)
blas_set_num_threads(1)
omp_set_num_threads(1)

setwd("/home/luc/Rpackages/MRAinla/")

rawDataFilesLocation <- "data-raw/"

### Preparing datasets...

# Naming convention: nnnnnnn.Ayyyyddd.h00v00.vvv.yyyydddhhmmss.
# nnnnnnn: Product name
# Ayyyyddd: Sampling date, year (yyyy), then day (ddd), between 1 and 365.
# h00v00: Identifies the grid tile (see https://lpdaac.usgs.gov/dataset_discovery/modis)
# vvv: Data version
# yyyydddhhmmss: Date data were processed, year, day, hour, minute, second.
dayOffset <- 121
dayRange <- 18:24
collectionDates <- paste("May", dayRange, "_2012", sep = "")
collectionDatesPOSIX <- as.POSIXct(paste("2012-05-", dayRange, sep = ""))

splitTemperaturesBySatellite <- lapply(c(Terra = "MOD11A1.A2012", Aqua = "MYD11A1.A2012"), function(searchString) {
  temperatureFiles <- list.files(path = rawDataFilesLocation, pattern = searchString, full.names = TRUE)
  subFiles <- sapply(paste("A2012", dayOffset + dayRange, sep = ""), grep, x = temperatureFiles, value = TRUE)
  temperatures <- lapply(subFiles, getSds)
  splitTemperatures <- split(temperatures, f = factor(substr(subFiles, start = 0, stop = gregexpr(pattern = ".h2", text = subFiles[[1]])[[1]] - 1)))
  names(splitTemperatures) <- collectionDates
  splitTemperatures
})

data(wrld_simpl)

# Mumbai polygon:

mumbaiPolygonEdges <- rbind(c(19.3, 72.76), c(19.3, 73.4), c(18.85, 73.4), c(18.85, 72.76), c(19.3, 72.76))
mumbaiPolygonEdges <- mumbaiPolygonEdges[ , 2:1]
mumbaiPolygon <- SpatialPolygons(Srl = list(Polygons(list(Polygon(coords = mumbaiPolygonEdges)), ID = "Mumbai")))
crs(mumbaiPolygon) <- crs(wrld_simpl)
mumbaiPolygonOtherCRS <- spTransform(mumbaiPolygon, CRSobj = crs(readGDAL(splitTemperaturesBySatellite$Aqua[[1]][[1]]$SDS4gdal[1], as.is = TRUE)))

indiaTemperaturesAndTimes <- lapply(seq_along(splitTemperaturesBySatellite$Aqua), function(var1) {
  aquaRasters <- funToCreateRaster(splitTemperaturesBySatellite$Aqua[[var1]], polygonBound = mumbaiPolygonOtherCRS)
  terraRasters <- funToCreateRaster(splitTemperaturesBySatellite$Terra[[var1]], polygonBound = mumbaiPolygonOtherCRS)
  if (sum(!is.na(values(aquaRasters$temperatureRaster))) >= sum(!is.na(values(terraRasters$temperatureRaster)))) {
    cat("Returning Aqua!\n")
    c(aquaRasters, satellite = "Aqua")
  } else {
    cat("Returning Terra!\n")
    c(terraRasters, satellite = "Terra")
  }
})

indiaTemperatures <- lapply(indiaTemperaturesAndTimes, function(x) x$temperatureRaster)
indiaTimes <- lapply(indiaTemperaturesAndTimes, function(x) x$hourRaster)
satellitePerDay <- sapply(indiaTemperaturesAndTimes, function(x) x$satellite)

# jpeg("outputFiles/temperaturesMay18dataExample.jpg", width = 1600, height = 1600)
# spplot(indiaTemperatures[[1]], scales = list(draw = TRUE), par.settings=list(fontsize=list(text=50)))
# dev.off()

# sum(sapply(indiaTemperatures, function(x) sum(!is.na(values(x)))))

landCoverFiles <- list.files(rawDataFilesLocation, pattern = "MCD*", full.names = TRUE)

produceLandCover <- function(landCoverFiles) { # Probably not absolutely necessary to subset land cover values, as the temperature values are already subsetted to be in India only.
  landCoverRasters <- lapply(landCoverFiles, function(filename) {
    landCoverSds <- getSds(filename)
    landCover <- raster(readGDAL(landCoverSds$SDS4gdal[2], as.is = TRUE)) # Based on land type classification 2: https://lpdaac.usgs.gov/products/mcd12q1v006/
    landCover
  })
  landCover <- do.call(raster::merge, landCoverRasters)
  # values(landCover) <- factor(values(landCover)) # Can't recall why I had included this line.
  smallerRaster <- crop(x = landCover, y = mumbaiPolygonOtherCRS) # This will only keep points in the box defined by mumbaiPolygonOtherCRS.
  # The next few lines keep only land cover values within an arbitrarily-shaped box.
  spObject <- rasterToPoints(smallerRaster, spatial = TRUE)
  indiaValuesIndex <- over(x = spObject, y = mumbaiPolygonOtherCRS)
  pointsInIndia <- subset(spObject, subset = !is.na(indiaValuesIndex))
  values(smallerRaster) <- rep(NA, ncell(smallerRaster))
  output <- rasterize(x = pointsInIndia, y = smallerRaster, field = "layer")
  output
}

landCover <- produceLandCover(landCoverFiles)

# Getting elevation data

elevationFiles <- list.files(path = rawDataFilesLocation, pattern = "*dem.tif", full.names = TRUE)
elevation <- lapply(elevationFiles, raster)

may21rasterAddedMissing <- indiaTemperatures[[4]]
may28raster <- indiaTemperatures[[length(indiaTemperatures) - 3]]
indicesForKnownRemovedValues <- which(is.na(values(may28raster)) & !is.na(values(may21rasterAddedMissing)))
values(may21rasterAddedMissing) <- replace(values(may21rasterAddedMissing), indicesForKnownRemovedValues, NA)

MODISdataMumbaiSinusoidalProj <- prepareDataForMRAinla(landCover = landCover, elevations = elevation, temperatures = indiaTemperatures, collectionDatesPOSIX = collectionDatesPOSIX, satelliteNamesVec = satellitePerDay) # Remember: covariates will have to be recentered prior to analyses.

MODISdataMumbai <- MODISdataMumbaiSinusoidalProj
MODISdataMumbai@sp <- sp::spTransform(x = MODISdataMumbai@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

modelRasterForPreds <- indiaTemperatures[[1]]

predRasterList <- lapply(indiaTemperatures, function(temperatureRaster) {
  rasterToFill <- modelRasterForPreds
  values(rasterToFill) <- rep(NA, ncell(rasterToFill))
  indicesForPredRaster <- which(!is.na(values(modelRasterForPreds)) & is.na(values(temperatureRaster)))
  if (length(indicesForPredRaster) == 0) return(NA)
  rasterToFill[indicesForPredRaster] <- 500
  rasterToFill
})
collectionIndices <- sapply(predRasterList, function(x) !identical(x, NA))
# sapply(predRasterList[collectionIndices], function(x) sum(!is.na(values(x))))
datesForPred <- collectionDatesPOSIX[collectionIndices]
satellitesforPred <- satellitePerDay[collectionIndices]
datesForPred <- collectionDatesPOSIX[collectionIndices]

MODISdataMumbaiTestSinusoidalProj <- prepareDataForMRAinla(landCover = landCover, elevations = elevation, temperatures = predRasterList[collectionIndices], collectionDatesPOSIX = datesForPred, satelliteNamesVec = satellitesforPred, completeDateVector = collectionDatesPOSIX) # Remember: covariates will have to be recentered prior to analyses.
MODISdataMumbaiTest <- MODISdataMumbaiTestSinusoidalProj
MODISdataMumbaiTest@sp <- sp::spTransform(x = MODISdataMumbaiTest@sp, CRSobj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

missingColumnNames <- setdiff(colnames(MODISdataMumbai@data), colnames(MODISdataMumbaiTest@data))
columnsToAdd <- lapply(missingColumnNames, function(columnName) rep(0, nrow(MODISdataMumbaiTest@data)))
# We need to add columns to the prediction dataset to have it match with the training dataset
names(columnsToAdd) <- missingColumnNames
columnsToAdd <- as.data.frame(columnsToAdd)
MODISdataMumbaiTest@data <- cbind(MODISdataMumbaiTest@data, columnsToAdd)
MODISdataMumbaiTest@data <- MODISdataMumbaiTest@data[ , colnames(MODISdataMumbai@data)]

# Water is the reference category.
MODISdataMumbai@data <- subset(MODISdataMumbai@data, select = -landCover0)
MODISdataMumbaiTest@data <- subset(MODISdataMumbaiTest@data, select = -landCover0)
MODISdataMumbaiTest@data <- subset(MODISdataMumbaiTest@data, select = -y)

MODISdataTraining <- MODISdataMumbai
MODISdataTest <- MODISdataMumbaiTest

usethis::use_data(MODISdataTraining, MODISdataTest, overwrite = TRUE)

