plot.INLAMRA <- function(x, filename = NULL, type = c("joint", "training", "predictions", "SD", "hyperpars"), polygonsToOverlay = NULL, control = plot.control()) {
  type <- type[[1]]

  .checkPlotInputConsistency(x, type)

  if (control$trim) {
    x$predMoments$Mean <- replace(x$predMoments$Mean, which(x$predMoments$Mean > max(x$data$spObject@data$y)), max(x$data$spObject@data$y))
    x$predMoments$Mean <- replace(x$predMoments$Mean, which(x$predMoments$Mean < min(x$data$spObject@data$y)), min(x$data$spObject@data$y))
  }

  if (control$plotRaster & (type != "hyperpars")) {
    .plotRaster(x, control, type = type)
  } else if (type[[1]] != "hyperpars") {
    .plotPoints(x, control, type = type)
  } else {
    .plotHyperparMarginal(output = x, filename = filename, device = control$graphicsDevice, width = control$width, height = control$height)
  }
}

.checkPlotInputConsistency <- function(INLAMRAoutput, type) {
  if (!(type %in% c("training", "predictions", "joint", "SD", "hyperpars"))) {
    stop("type argument must be one of training, predictions, joint, SD.")
  }
  if (!(type == "hyperpars") & is.null(INLAMRAoutput$predData)) {
    stop("Requested to get prediction graphs. If no prediction data were provided, only marginal hyperparameter posterior graphs are available. Please set type to 'hyperpars'.")
  }
  invisible(0)
}

.plotRaster <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type) {
  testDays <- sort(unique(INLAMRAoutput$predData$time))
  control <- .SetRasterSizes(control, INLAMRAoutput, type = type)
  rasterListPerDay <- lapply(testDays, .rasterizeTrainingTestJointSD, INLAMRAoutput = INLAMRAoutput, control = control)

  if (!is.null(filename)) {
    control$graphicsEngine(filename, width = control$width, height = control$height)
  }

  rastersToPlot <- lapply(rasterListPerDay, function(x) x[[type]])
  rasterRanges <- sapply(rastersToPlot, function(x) range(raster::values(x)))
  colorRange <- c(min(rasterRanges[1, ]), max(rasterRanges[2, ]))

  if (control$matchColours) {
    rastersForRange <- lapply(rasterListPerDay, function(x) x[["joint"]])
    rasterRanges <- sapply(rastersForRange, function(x) range(raster::values(x)))
    colorRange <- c(min(rasterRanges[1, ]), max(rasterRanges[2, ]))
  }

  numPlotsPerLine <- ceiling(sqrt(length(rastersToPlot)))
  par(mfrow = c(numPlotsPerLine, numPlotsPerLine), mai = rep(1.5, 4))
  for (i in seq_along(rastersToPlot)) {
    raster::plot(rastersToPlot[[i]], zlim = colorRange)
    raster::plot(rastersToPlot[[i]], legend.only = TRUE, legend.width = 5, axis.args = list(cex.axis = 4), zlim = colorRange)
    if (!is.null(polygonsToOverlay)) {
      raster::plot(polygonsToOverlay, add = TRUE)
      mapmisc::scaleBar(crs = raster::crs(polygonsToOverlay), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
    }
  }
  if (!is.null(filename)) {
    dev.off()
  }
  invisible(0)
}

.SetRasterSizes <- function(control, INLAMRAoutput, type) {
  uniqueTimeValues <- sort(unique(c(INLAMRAoutput$data$time, INLAMRAoutput$predData$time)))

  findNrowsNcolsByTimeIndex <- function(timeValue, coordIndex) {
    coordsForTraining <- INLAMRAoutput$data$spObject@coords[INLAMRAoutput$data$time == timeValue, ]
    coordsForPredSet <- INLAMRAoutput$predData$spObject@coords[INLAMRAoutput$predData$time == timeValue, ]
    if (type == "joint") {
      combinedCoords <- rbind(coordsForTraining, coordsForPredSet)
    } else if (type == "training") {
      combinedCoords <- coordsForTraining
    } else {
      combinedCoords <- coordsForPredSet
    }
    orderedCombinedCoords <- combinedCoords[order(combinedCoords[ , coordIndex]), ]
    diffCoord <- diff(orderedCombinedCoords[, coordIndex])
    sdDiffCoord <- sd(diffCoord)
    breakPositions <- which(abs(diffCoord) > 3*sdDiffCoord)
    dimSizes <- c(breakPositions[[1]], diff(breakPositions))
    max(dimSizes)
  }
  proposedNrowsNcolsByTime <- sapply(uniqueTimeValues, FUN = function(timeValue) {
    sapply(c(2,1),  FUN = findNrowsNcolsByTimeIndex, timeValue = timeValue)
  })

  control$rasterNrows <- max(proposedNrowsNcolsByTime[1, ])
  control$rasterNcols <- max(proposedNrowsNcolsByTime[2, ])
  cat("Trying to infer the ideal number of cells in the raster by rounding spatial coordinates at", control$numDigitRound,"digits to eliminate the spatial jittering created by INLAMRA. Obtained", control$rasterNrows, "rows and", control$rasterNcols, "columns. If the raster looks bad in the end, set these values manually with plot.control(). If your data are not gridded, set plotRaster to FALSE in plot.control() instead.", sep = " ")
  control
}

# Only applied to time points where predictions were required.

.rasterizeTrainingTestJointSD <- function(timePoint, INLAMRAoutput, landRaster, control) {
  combinedCoordMat <- rbind(INLAMRAoutput$data$spObject@coords, INLAMRAoutput$predData$spObject@coords)
  landRaster <- raster::raster(nrows = control$rasterNrows, ncols = control$rasterNcols, xmn = min(combinedCoordMat[ , 1]), xmx = max(combinedCoordMat[ , 1]), ymn = min(combinedCoordMat[ , 2]), ymx = max(combinedCoordMat[ , 2]), crs = sp::proj4string(INLAMRAoutput$data$spObject))

  trainingDataIndices <- INLAMRAoutput$data$time == timePoint
  testDataIndices <- INLAMRAoutput$predData$time == timePoint

  rasterList <- list()
  rasterList$SD <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$SD[testDataIndices])
  rasterList$predictions <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$Mean[testDataIndices])
  if (any(trainingDataIndices)) {
    rasterList$training <- raster::rasterize(x = INLAMRAoutput$data$spObject@coords[trainingDataIndices, ], y = landRaster, field = INLAMRAoutput$data$spObject@data[trainingDataIndices , "y"])

    jointCoordinates <- unname(rbind(INLAMRAoutput$predData$spObject@coords[testDataIndices, ], INLAMRAoutput$data$spObject@coords[trainingDataIndices, ]))
    dataObject <- data.frame(y = unname(c(INLAMRAoutput$predMoments$Mean[testDataIndices], INLAMRAoutput$data$spObject@data[trainingDataIndices, "y"])))
    pointsDataFrame <- sp::SpatialPointsDataFrame(coords = jointCoordinates, data = dataObject)
    rasterList$joint <- raster::rasterize(x = pointsDataFrame, y = landRaster, field = "y")
  } else {
    rasterList$training <- landRaster
    rasterList$joint <- rasterList$predictions
  }
  rasterList
}

.plotPoints <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type) {
  pointsPerDayList <- .getPointsPerDay(INLAMRAoutput = INLAMRAoutput, control = control)
  pointsToPlotList <- lapply(pointsPerDayList, function(x) x[[type]])

  valuesRanges <- sapply(pointsToPlotList, function(x) range(x@data$y))
  colorRange <- c(min(valuesRanges[1, ]), max(valuesRanges[2, ]))

  if (control$matchColours) {
    pointsForRange <- lapply(pointsToPlotList, function(x) x[["joint"]])
    pointColourRanges <- sapply(pointsForRange, function(x) range(x@data$y))
    colorRange <- c(min(pointColourRanges[1, ]), max(pointColourRanges[2, ]))
  }

  numPlotsPerLine <- ceiling(sqrt(length(pointsToPlotList)))
  par(mfrow = c(numPlotsPerLine, numPlotsPerLine), mai = rep(1.5, 4))

  for (i in seq_along(pointsToPlotList)) {
    sp::plot(pointsToPlotList[[i]], zlim = colorRange)
    sp::plot(pointsToPlotList[[i]], legend.only = TRUE, legend.width = 5, axis.args = list(cex.axis = 4), zlim = colorRange)
    if (!is.null(polygonsToOverlay)) {
      sp::plot(polygonsToOverlay, add = TRUE)
      mapmisc::scaleBar(crs = raster::crs(polygonsToOverlay), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
    }
  }
  if (!is.null(filename)) {
    dev.off()
  }
  invisible(0)
}

.getPointsPerDay <- function(INLAMRAoutput, control) {
  testDays <- sort(unique(INLAMRAoutput$predData$time))
  getPoints <- function(day) {
    trainingPointsIndices <- INLAMRAoutput$time == day
    trainingPoints <- NULL
    if (any(trainingPointsIndices)) {
      trainingPoints <- INLAMRAoutput$data$spObject[trainingPointsIndices]
    }
    predPointsIndices <- INLAMRAoutput$predData == day
    predPoints <- INLAMRAoutput$predData$spObject[predPointsIndices]
    predPoints@data <- cbind(data.frame(y = INLAMRAoutput$predMoments$Mean), predPoints@data)
    jointPoints <- c(trainingPoints, predPoints)
    list(training = trainingPoints, predictions = predPoints, joint = jointPoints)
  }
  lapply(testDays, getPoints)
}

plot.control <- function(trim = FALSE, fontScaling = 1, width = 1000, height = width, plotRaster = TRUE, rasterNrows = NULL, rasterNcols = NULL, numDigitRound = 5, graphicsEngine = tiff, matchColours = FALSE) {
  list(trim = trim, fontScaling = fontScaling, width = width, height = height, plotRaster = plotRaster, rasterNrows = rasterNrows, rasterNcols = rasterNcols, numDigitRound = numDigitRound, graphicsEngine = graphicsEngine, matchColours = matchColours)
}

.plotHyperparMarginal <- function(output, numValues = 50, device = jpeg, filename, ...) {
  plotValuesList <- function(parameterName) {
    hyperparSkewness <- output$hyperMarginalMoments[parameterName, "Skewness"]
    hyperparSD <- output$hyperMarginalMoments[parameterName, "StdDev"]
    hyperparMean <- output$hyperMarginalMoments[parameterName, "Mean"]
    if (!(is.na(hyperparSD) | (hyperparSD > 0))) {
      credIntBounds <- .ComputeCredIntervalSkewNorm(c(0.025, 0.975), meanValue = hyperparMean, sdValue = hyperparSD, skewnessValue = hyperparSkewness)
    } else {
      stop("Error: selected hyperparameter was fixed when INLA-MRA was used.")
    }
    xValues <- seq(from = credIntBounds$bounds[1], to = credIntBounds$bounds[2], length.out = numValues)
    yValues <- sn::dsn(x = xValues, xi = credIntBounds$xi, omega = credIntBounds$omega, alpha = credIntBounds$alpha)
    plotFrame <- data.frame(x = xValues, Value = yValues)
    colnames(plotFrame)[[1]] <- parameterName
    plotFrame
  }
  plotFrames <- lapply(rownames(output$hyperMarginalMoments), plotValuesList)
  par(mfrow = c(ceiling(sqrt(length(plotFrames))), ceiling(sqrt(length(plotFrames)))))
  device(file = filename, ...)
  for (i in seq_along(plotFrames)) {
    plot(x = plotFrames[[i]][ , 1], y = plotFrames[[i]]$Value, xlab = colnames(plotFrames[[i]])[[1]], ylab = "Value", type = "l")
  }
  dev.off()
}
