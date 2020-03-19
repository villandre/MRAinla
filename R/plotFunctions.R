#' Plots the output from INLAMRA
#'
#' Plots the output from INLAMRA. It prints the predictions by default. Note that this is only possible if the data were saved in the output (which is the case by default).
#'
#' @param x INLAMRA object. Created after calling `INLAMRA`
#' @param filename Character. File in which to save the output.
#' @param type Character. The default, "joint" creates a grid of plots, one plot per time value, showing the training data, and the training data incremented with the predictions. ""SD" creates a grid of plots showing standard deviations for predictions. "trainingData" and "predictions" produce a similar grid with only the training data and predictions, respectively. "pars" produces plots of the marginal posterior densities for the parameters and hyperparameters.
#' @param polygonsToOverlay Optional SpatialPolygons object. Those polygons, e.g. regional boundaries, are added to each prediction map.
#' @param ... Arguments for the graphics engine, e.g. jpeg. 'width' and 'height' should be especially helpful.
#' @param control Output of `plot.control`. See ?plot.control.
#'
#' @details The function produces a grid of rasters by default, cf. plot.control, one per distinct time value. If your data are not spaced regularly, you can plot the result as points instead by setting `control = plot.control(plotRaster = FALSE)`.
#'
#'
#' @return A list of the plotted raster/SpatialPoints objects if type is anything but "pars"; a list of matching x/y values for the plotted (hyper)parameter posteriors otherwise.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export
#'
plot.INLAMRA <- function(x, filename = NULL, type = c("joint", "training", "predictions", "SD", "marginals"), polygonsToOverlay = NULL, control = plot.control(), ...) {
  type <- type[[1]]

  .checkPlotInputConsistency(x, type)

  if (control$trim) {
    x$predMoments$Mean <- replace(x$predMoments$Mean, which(x$predMoments$Mean > max(x$data$spObject@data$y)), max(x$data$spObject@data$y))
    x$predMoments$Mean <- replace(x$predMoments$Mean, which(x$predMoments$Mean < min(x$data$spObject@data$y)), min(x$data$spObject@data$y))
  }

  if (control$plotRaster & (type != "marginals")) {
    plottedObjects <- .plotRaster(x, control, filename = filename, type = type, polygonsToOverlay = polygonsToOverlay, ...)
  } else if (type[[1]] != "marginals") {
    plottedObjects <- .plotPoints(x, control, filename = filename, type = type, polygonsToOverlay = polygonsToOverlay, ...)
  } else {
    plottedObjects <- .plotMarginals(output = x, filename = filename, device = control$graphicsEngine, ...)
  }
  plottedObjects
}

.checkPlotInputConsistency <- function(INLAMRAoutput, type) {
  if (!(type %in% c("training", "predictions", "joint", "SD", "marginals"))) {
    stop("type argument must be one of training, predictions, joint, SD.")
  }
  if (!(type == "marginals") & is.null(INLAMRAoutput$predData)) {
    stop("Requested to get prediction graphs. If no prediction data were provided, only marginal hyperparameter posterior graphs are available. Please set type to 'hyperpars'.")
  }
  invisible(0)
}

.plotRaster <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type, ...) {
  testDays <- sort(unique(INLAMRAoutput$predData$time))
  control <- .SetRasterSizes(control, INLAMRAoutput, type = type)
  rasterListPerDay <- lapply(testDays, .rasterizeTrainingTestJointSD, INLAMRAoutput = INLAMRAoutput, control = control)

  rastersToPlot <- lapply(rasterListPerDay, function(x) x[[type]])
  rasterRanges <- sapply(rastersToPlot, function(x) range(raster::values(x)))
  colorRange <- c(min(rasterRanges[1, ]), max(rasterRanges[2, ]))

  if (control$matchColours) {
    rastersForRange <- lapply(rasterListPerDay, function(x) x[["joint"]])
    rasterRanges <- sapply(rastersForRange, function(x) range(raster::values(x)))
    colorRange <- c(min(rasterRanges[1, ]), max(rasterRanges[2, ]))
  }

  numPlotsPerLine <- ceiling(sqrt(length(rastersToPlot)))
  if (!is.null(filename)) {
    control$graphicsEngine(filename = filename, ...)
  }
  layout(matrix(1:(numPlotsPerLine^2), nrow = numPlotsPerLine, ncol = numPlotsPerLine, byrow = TRUE))
  for (i in seq_along(rastersToPlot)) {
    raster::plot(rastersToPlot[[i]], zlim = colorRange, interpolate = FALSE, useRaster = FALSE, legend = FALSE)
    raster::plot(rastersToPlot[[i]], legend.only = TRUE, legend.width = 4, axis.args = list(cex.axis = 3), zlim = colorRange)
    if (!is.na(raster::crs(rastersToPlot[[i]]))) {
      mapmisc::scaleBar(crs = raster::crs(rastersToPlot[[i]]), pos = "topleft", cex = 2, pt.cex = 1.2, title.cex = 1.2)
    }
    if (!is.null(polygonsToOverlay)) {
      raster::plot(polygonsToOverlay, add = TRUE)
    }
  }
  if (!is.null(filename)) {
    dev.off()
  }
  rastersToPlot
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

  control$rasterNrows <- max(proposedNrowsNcolsByTime[1, ])  # The multiplier is there to remove white lines in the raster when the number of rows is not estimated perfectly.
  control$rasterNcols <- max(proposedNrowsNcolsByTime[2, ])  # The multiplier is there to remove white lines in the raster when the number of rows is not estimated perfectly.
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

.plotPoints <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type, ...) {
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

  if (!is.null(filename)) {
    control$graphicsEngine(filename = filename, ...)
  }
  layout(matrix(1:(numPlotsPerLine^2), nrow = numPlotsPerLine, ncol = numPlotsPerLine, byrow = TRUE))
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
  pointsToPlotList
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

plot.control <- function(trim = FALSE, fontScaling = 1, plotRaster = TRUE, rasterNrows = NULL, rasterNcols = NULL, numDigitRound = 5, graphicsEngine = jpeg, matchColours = FALSE) {
  list(trim = trim, fontScaling = fontScaling, plotRaster = plotRaster, rasterNrows = rasterNrows, rasterNcols = rasterNcols, numDigitRound = numDigitRound, graphicsEngine = graphicsEngine, matchColours = matchColours)
}

.plotMarginals <- function(output, numValues = 50, device = jpeg, filename, ...) {
  plotValuesList <- function(parameterName) {
    hyperparSkewness <- output$hyperMarginalMoments[parameterName, "Skewness"]
    hyperparSD <- output$hyperMarginalMoments[parameterName, "StdDev"]
    hyperparMean <- output$hyperMarginalMoments[parameterName, "Mean"]
    credIntBounds <- .ComputeCredIntervalSkewNorm(c(0.025, 0.975), meanValue = hyperparMean, sdValue = hyperparSD, skewnessValue = hyperparSkewness)
    xValues <- seq(from = credIntBounds$bounds[1], to = credIntBounds$bounds[2], length.out = numValues)
    yValues <- sn::dsn(x = xValues, xi = credIntBounds$xi, omega = credIntBounds$omega, alpha = credIntBounds$alpha)
    data.frame(x = xValues, y = yValues)
  }
  hyperPlotFrames <- lapply(rownames(output$hyperMarginalMoments), plotValuesList)
  names(hyperPlotFrames) <- rownames(output$hyperMarginalMoments)
  plotFrames <- c(hyperPlotFrames, output$FEmargDistValues)
  device(file = filename, ...)
  layout(matrix(1:ceiling(sqrt(length(plotFrames)))^2, nrow = ceiling(sqrt(length(plotFrames))), ncol = ceiling(sqrt(length(plotFrames))), byrow = TRUE))
  for (i in seq_along(plotFrames)) {
    plot(x = plotFrames[[i]]$x, y = plotFrames[[i]]$y, xlab = names(plotFrames)[[i]], ylab = "Value", type = "l", cex = 2)
  }
  dev.off()
}
