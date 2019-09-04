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

plotOutput <- function(inlaMRAoutput, trainingData, testData, realTestValues = NULL, filename = NULL, graphicsEngine = tiff, plotWhat = c("joint", "training", "SD", "fittedVSrealNoSp"), control = list()) {
  if (!("width" %in% names(control))) {
    control$width <- control$height <- 1600
  }

  if (!("rasterNrow" %in% names(control))) {
    control$rasterNrow <- 50
    control$rasterNcol <- 50
  }
  lonLatMinMax <- lapply(1:2, function(colIndex) {
    sapply(list(min, max), function(summaryFunction) {
      summaryFunction(testData@sp@coords[ , colIndex], trainingData@sp@coords[ , colIndex])
    })
  })

  landRaster <- raster::raster(nrows = control$rasterNrow, ncols = control$rasterNcol, xmn = lonLatMinMax[[1]][[1]], xmx = lonLatMinMax[[1]][[2]], ymn = lonLatMinMax[[2]][[1]], ymx = lonLatMinMax[[2]][[2]])
  uniqueTimeValues <- unique(c(time(trainingData), time(testData)))

  rasterizeTrainingAndJoint <- function(timePoint) {
    rasterList <- list()
    trainingDataIndices <- time(trainingData) == timePoint
    testDataIndices <- time(testData) == timePoint
    rasterList$landRasterTraining <- raster::rasterize(x = trainingData@sp@coords[trainingDataIndices, ], y = landRaster, field = trainingData@data[trainingDataIndices , 1])
    landRasterJointSD <- NULL
    landRasterTest <- NULL
    if (any(testDataIndices)) {
      landRasterNoSpFillMat <- matrix(rep(NA, control$rasterNrow * control$rasterNcol), control$rasterNrow, control$rasterNcol)
      rasterList$landRasterTestNoSp <- raster::raster(x = replace(landRasterNoSpFillMat, 1:sum(testDataIndices), realTestValues[testDataIndices]))
      rasterList$landRasterFittedNoSp <- raster::raster(x = replace(landRasterNoSpFillMat, 1:sum(testDataIndices), inlaMRAoutput$predictionMoments$predictMeans[testDataIndices]))
      rasterList$landRasterJointSD <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictSDs[testDataIndices])
      rasterList$landRasterFitted <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictMeans[testDataIndices])
    }
    if (any(testDataIndices) & any(trainingDataIndices)) {
      jointCoordinates <- unname(rbind(testData@sp@coords[testDataIndices, ], trainingData@sp@coords[trainingDataIndices, ]))
      dataObject <- data.frame(y = unname(c(inlaMRAoutput$predictionMoments$predictMeans[testDataIndices], trainingData@data[trainingDataIndices, 1])))
      pointsDataFrame <- sp::SpatialPointsDataFrame(coords = jointCoordinates, data = dataObject)
      rasterList$landRasterJoint <- raster::rasterize(x = pointsDataFrame, y = landRaster, field = "y")
      rasterList$landRasterTest <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictMeans[testDataIndices])
    } else if (any(trainingDataIndices)) {
      rasterList$landRasterJoint <- rasterList$landRasterTraining
    } else if (any(testDataIndices)) {
      rasterList$landRasterJoint <- raster::rasterize(x = testData@sp@coords[testDataIndices, ], y = landRaster, field = inlaMRAoutput$predictionMoments$predictMeans[testDataIndices])
    } else {
      rasterList$landRasterJoint <- NULL
    }
    list(training = rasterList$landRasterTraining, joint = rasterList$landRasterJoint, fitted = rasterList$landRasterFitted, SD = rasterList$landRasterJointSD, testNoSp = rasterList$landRasterTestNoSp, fittedNoSp = rasterList$landRasterFittedNoSp)
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
    graphicsEngine(filename, width = control$width, height = control$height)
  }
  if (plotWhat == "SD") {
    stackedRasters <- raster::stack(stackedRastersList$fitted, stackedRastersList$SD)
  } else if (plotWhat == "fittedVsRealNoSp") {
    stackedRasters <- raster::stack(stackedRastersList$testNoSp, stackedRastersList$fittedNoSp)
  } else if (plotWhat == "training") {
    stackedRasters <- raster::stack(stackedRastersList$training)
  } else if (plotWhat == "joint") {
    stackedRasters <- raster::stack(stackedRastersList$training, stackedRastersList$joint)
  } else {
    stop("Unrecognised plot type requested: please select one of: joint, training, SD, fittedVSrealNoSp.")
  }
  rangeForScale <- range(raster::values(stackedRasters), na.rm = TRUE)
  plot(stackedRasters, interpolate = FALSE, col = rev( rainbow( 20, start = 0, end = 1) ), breaks = seq(floor(rangeForScale[[1]]), ceiling(rangeForScale[[2]]), length.out = 19))

  if (!is.null(filename)) {
    dev.off()
  }
  invisible(0)
}
