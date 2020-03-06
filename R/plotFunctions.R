plot.INLAMRA <- function(INLAMRAoutput, filename = NULL, type = c("joint", "training", "SD"), polygonsToOverlay = NULL, control = plot.control()) {
  .checkPlotInputConsistency(INLAMRAoutput)

  if (control$trim) {
    INLAMRAoutput$predMoments$Mean <- replace(INLAMRAoutput$predMoments$Mean, which(INLAMRAoutput$predMoments$Mean > max(INLAMRAoutput$data$spObject@data$y)), max(INLAMRAoutput$data$spObject@data$y))
    INLAMRAoutput$predMoments$Mean <- replace(INLAMRAoutput$predMoments$Mean, which(INLAMRAoutput$predMoments$Mean < min(INLAMRAoutput$data$spObject@data$y)), min(INLAMRAoutput$data$spObject@data$y))
  }

  if (control$plotRaster) {
    .plotRaster(INLAMRAoutput, control, type = type)
  } else {
    .plotPoints(INLAMRAoutput, control, type = type)
  }
}

.plotRaster <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type) {
  testDays <- sort(unique(INLAMRAoutput$predData$time))
  control <- .SetRasterSizes(control, INLAMRAoutput)
  uniqueTimeValues <- testDays
  dailyRasters <- lapply(uniqueTimeValues, .rasterizeTrainingAndJoint, INLAMRAoutput = INLAMRAoutput, control = control)

  rasterNames <- names(dailyRasters[[1]])
  stackedRastersList <- lapply(rasterNames, FUN = .funToGetStackedRaster, dailyRasters = dailyRasters, uniqueTimeValues = uniqueTimeValues)
  names(stackedRastersList) <- rasterNames

  if (!is.null(filename)) {
    control$graphicsEngine(filename, width = control$width, height = control$height)
  }
  if (type == "SD") {
    stackedRasters <- stackedRastersList$SD
  } else if (type == "training") {
    stackedRasters <- stackedRastersList$training
  } else if (type == "joint") {
    stackedRasters <- c(stackedRastersList$training, stackedRastersList$joint)
  } else {
    stop("Unrecognised plot type requested: please select one of: joint, training, SD, fittedVSrealNoSp.")
  }
  par(mfrow = c(2, 3), mai = rep(1.5, 4))
  for (i in seq_along(stackedRasters)) {
    raster::plot(stackedRasters[[i]], legend = FALSE, axes = FALSE)
    if (!is.null(polygonsToOverlay)) {
      raster::plot(polygonsToOverlay, add = TRUE)
      mapmisc::scaleBar(crs = raster::crs(polygonsToOverlay), pos = "topleft", cex = 3, pt.cex = 2.2, title.cex = 3.5)
    }
    raster::plot(stackedRasters[[i]], legend.only = TRUE, legend.width = 5, axis.args = list(cex.axis = 4))
  }
  # plot(stackedRasters, interpolate = FALSE, col = rev( rainbow( 20, start = 0, end = 1) ), breaks = seq(floor(rangeForScale[[1]]), ceiling(rangeForScale[[2]]), length.out = 19), cex = control$fontScaling)

  if (!is.null(filename)) {
    dev.off()
  }
  invisible(0)
}

.SetRasterSizes <- function(control, INLAMRAoutput) {
  combinedCoords <- rbind(INLAMRAoutput$data$spObject@coords, INLAMRAoutput$predData$spObject@coords)
  # We try to remove the effects of the jittering.
  control$rasterNrows <- length(unique(round(combinedCoords[ , 1], digits = control$numDigitRound)))
  control$rasterNcols <- length(unique(round(combinedCoords[ , 2], digits = control$numDigitRound)))
  cat("Trying to infer the ideal number of cells in the raster by rounding spatial coordinates at", control$numDigitRound,"digits to eliminate the spatial jittering created by INLAMRA. Obtained", control$rasterNrows, "rows and", control$rasterNcols, "columns. If the raster looks bad in the end, set these values manually with plot.control(). If your data are not gridded, set plotRaster to FALSE in plot.control() instead.", sep = " ")
  control
}

.rasterizeTrainingAndJoint <- function(timePoint, INLAMRAoutput, landRaster, control) {
  combinedCoordMat <- rbind(INLAMRAoutput$data$spObject@coords, INLAMRAoutput$predData$spObject@coords)
  landRaster <- raster::raster(nrows = control$rasterNrows, ncols = control$rasterNcols, xmn = min(combinedCoordMat[ , 1]), xmx = max(combinedCoordMat[ , 1]), ymn = min(combinedCoordMat[ , 2]), ymx = max(combinedCoordMat[ , 2]), crs = sp::proj4string(INLAMRAoutput$data$spObject))

  rasterList <- list()
  trainingDataIndices <- INLAMRAoutput$data$time == timePoint
  testDataIndices <- INLAMRAoutput$predData$time == timePoint
  if (any(trainingDataIndices)) {
    rasterList$landRasterTraining <- raster::rasterize(x = INLAMRAoutput$data$spObject@coords[trainingDataIndices, ], y = landRaster, field = INLAMRAoutput$data$spObject@data[trainingDataIndices , "y"])
  } else {
    rasterList$landRasterTraining <- landRaster
  }

  rasterList$landRasterJointSD <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$SD[testDataIndices])
  rasterList$landRasterFitted <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$Mean[testDataIndices])

  if (any(trainingDataIndices)) {
    jointCoordinates <- unname(rbind(INLAMRAoutput$predData$spObject@coords[testDataIndices, ], INLAMRAoutput$data$spObject@coords[trainingDataIndices, ]))
    dataObject <- data.frame(y = unname(c(INLAMRAoutput$predMoments$Mean[testDataIndices], INLAMRAoutput$data$spObject@data[trainingDataIndices, "y"])))
    pointsDataFrame <- sp::SpatialPointsDataFrame(coords = jointCoordinates, data = dataObject)
    rasterList$landRasterJoint <- raster::rasterize(x = pointsDataFrame, y = landRaster, field = "y")
    rasterList$landRasterTest <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$Mean[testDataIndices])
  } else if (any(trainingDataIndices)) {
    rasterList$landRasterJoint <- rasterList$landRasterTraining
  } else {
    rasterList$landRasterJoint <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$Mean[testDataIndices])
  }
  output <- list(training = rasterList$landRasterTraining, joint = rasterList$landRasterJoint, fitted = rasterList$landRasterFitted, SD = rasterList$landRasterJointSD)

  output
}

.funToGetStackedRaster <- function(dataName, dailyRasters, uniqueTimeValues) {
  nullAndEmptyPos <- sapply(dailyRasters, function(x) {
    if (is.null(x[[dataName]])) {
      return(TRUE)
    } else {
      if (all(is.na(raster::values(x[[dataName]])))) return(TRUE)
    }
    FALSE
  })
  stackedRasters <- lapply(dailyRasters, function(x) x[[dataName]])[!nullAndEmptyPos]
  names(stackedRasters) <- paste(dataName, ":", as.character(uniqueTimeValues[!nullAndEmptyPos]), sep = "")
  stackedRasters
}

.plotPoints <- function(INLAMRAoutput, control) {

}

plot.control <- function(trim = FALSE, fontScaling = 1, width = 1000, height = width, plotRaster = TRUE, rasterNrows = NULL, rasterNcols = NULL, numDigitRound = 5, graphicsEngine = tiff) {
  list(trim = trim, fontScaling = fontScaling, width = width, height = height, plotRaster = plotRaster, rasterNrows = rasterNrows, rasterNcols = rasterNcols, numDigitRound = numDigitRound, graphicsEngine = graphicsEngine)
}

.plotHyperparMarginal <- function(output, parameterName, numValues, device, filename, ...) {
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
  device(file = filename, ...)
  plot(x = xValues, y = yValues, xlab = parameterName, ylab = "Value", type = "l")
  dev.off()
}
