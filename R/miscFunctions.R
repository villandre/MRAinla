SimulateSpacetimeData <- function(numObsPerTimeSlice = 225, covFunction, lonRange, latRange, timeValuesInPOSIXct, covariateGenerationFctList, errorSD, distFun = dist, FEvalues) {
  numSlotsPerRow <- ceiling(sqrt(numObsPerTimeSlice))
  slotCoordinates <- sapply(list(longitude = lonRange, latitude = latRange), FUN = function(x) {
    width <- abs(diff(x)/numSlotsPerRow)
    seq(from = min(x) + width/2, to = max(x) - width/2, length.out = numSlotsPerRow)
  }, USE.NAMES = TRUE)

  allSpaceCoordinates <- as.data.frame(expand.grid(as.data.frame(slotCoordinates)))
  numToRemove <- nrow(allSpaceCoordinates) - numObsPerTimeSlice

  if (numToRemove > 0) {
    obsToRemove <- (nrow(allSpaceCoordinates) - numToRemove + 1):nrow(allSpaceCoordinates)
    allSpaceCoordinates <- allSpaceCoordinates[-obsToRemove, ]
  }

  coordinates <- allSpaceCoordinates[rep(1:nrow(allSpaceCoordinates), length(timeValuesInPOSIXct)), ]
  coordinates$time <- rep(timeValuesInPOSIXct, each = numObsPerTimeSlice)

  covariateMatrix <- cbind(1, do.call("cbind", lapply(covariateGenerationFctList, function(x) {
    covariateValues <- x(coordinates[, c("longitude", "latitude")], coordinates[, "time"])
    if (is.null(dim(covariateValues))) {
      dim(covariateValues) <- c(length(covariateValues), 1)
    }
    covariateValues
  })))

  spatialDistMatrix <- dist(coordinates[, c("longitude", "latitude")])
  timeDistMatrix <- dist(coordinates[, "time"])/(3600*24)
  covarianceMat <- covFunction(spatialDistMatrix, timeDistMatrix)
  meanVector <- drop(covariateMatrix %*% FEvalues)
  fieldValues <- drop(mvtnorm::rmvnorm(n = 1, mean = meanVector, sigma = covarianceMat)) + rnorm(n = length(meanVector), mean = 0, sd = errorSD)
  dataForObject <- cbind(y = fieldValues, as.data.frame(covariateMatrix[, -1]))
  colnames(dataForObject) <- c("y", paste("Covariate", 1:(length(FEvalues) - 1), sep = ""))
  spacetimeObj <- spacetime::STIDF(sp = sp::SpatialPoints(coordinates[, c("longitude", "latitude")]), data = dataForObject, time = coordinates$time)
  spacetimeObj
}

plotOutput <- function(inlaMRAoutput, trainingData, testData, rasterNrow, rasterNcol, filename = NULL, graphicsEngine = tiff, plotWhat = c("joint", "training", "SD")) {
  lonLatMinMax <- lapply(1:2, function(colIndex) {
    sapply(list(min, max), function(summaryFunction) {
      summaryFunction(testData@sp@coords[ , colIndex], trainingData@sp@coords[ , colIndex])
    })
  })

  landRaster <- raster::raster(nrows = rasterNrow, ncols = rasterNcol, xmn = lonLatMinMax[[1]][[1]], xmx = lonLatMinMax[[1]][[2]], ymn = lonLatMinMax[[2]][[1]], ymx = lonLatMinMax[[2]][[2]])
  uniqueTimeValues <- unique(c(time(trainingData), time(testData)))

  rasterizeTrainingAndJoint <- function(timePoint) {
    trainingDataIndices <- time(trainingData) == timePoint
    testDataIndices <- time(testData) == timePoint
    landRasterTraining <- raster::rasterize(x = trainingData@sp@coords[trainingDataIndices, ], y = landRaster, field = trainingData@data[trainingDataIndices , 1])
    landRasterJointSD <- NULL
    landRasterTest <- NULL
    if (any(testDataIndices)) {
      landRasterJointSD <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictSDs)
    }
    if (any(testDataIndices) & any(trainingDataIndices)) {
      jointCoordinates <- unname(rbind(testData@sp@coords[testDataIndices, ], trainingData@sp@coords[trainingDataIndices, ]))
      dataObject <- data.frame(y = unname(c(inlaMRAoutput$predictionMoments$predictMeans[testDataIndices], trainingData@data[trainingDataIndices, 1])))
      pointsDataFrame <- sp::SpatialPointsDataFrame(coords = jointCoordinates, data = dataObject)
      landRasterJoint <- raster::rasterize(x = pointsDataFrame, y = landRaster, field = "y")
      landRasterTest <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictMeans[testDataIndices])
    } else if (any(trainingDataIndices)) {
      landRasterJoint <- landRasterTraining
    } else if (any(testDataIndices)) {
      landRasterJoint <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictMeans[testDataIndices])
    } else {
      landRasterJoint <- NULL
    }
    list(training = landRasterTraining, joint = landRasterJoint, SD = landRasterJointSD)
  }
  dailyRasters <- lapply(uniqueTimeValues, rasterizeTrainingAndJoint)

  rasterNames <- names(dailyRasters[[1]])
  stackedRastersList <- lapply(rasterNames, function(dataName) {
    nullPos <- sapply(dailyRasters, function(x) is.null(x[[dataName]]))
    stackedRasters <- raster::stack(lapply(dailyRasters, function(x) x[[dataName]])[!nullPos])
    names(stackedRasters) <- paste(dataName, ":", as.character(uniqueTimeValues[!nullPos]), sep = "")
    stackedRasters
  })
  names(stackedRastersList) <- rasterNames

  if (!is.null(filename)) {
    graphicsEngine(filename, width = 1600, height = 1600)
  }
  stackedRasters <- stack(stackedRastersList$training, stackedRastersList$joint, stackedRastersList$test)
  rasterRanges <- sapply(as.list(stackedRasters), FUN = function(aRasterLayer) range(raster::values(aRasterLayer), na.rm = TRUE))
  rangeForScale <- range(rasterRanges)
  plot(stackedRasters, interpolate = TRUE, col = rev( rainbow( 20, start = 0, end = 1) ), breaks = seq(floor(rangeForScale[[1]]), ceiling(rangeForScale[[2]]), length.out = 19))

  if (!is.null(filename)) {
    dev.off()
  }
  invisible(0)
}
