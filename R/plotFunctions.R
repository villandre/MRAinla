#' Plots the output from INLAMRA
#'
#' Plots the predictions or the marginal posteriors produced by INLAMRA. The predictions can only be plotted if the data are found in the INLAMRA object (it will be the case unless the user explicitly disabled it by setting control parameter saveData to FALSE).
#'
#' @param x INLAMRA object,
#' @param filename Character. File in which to save the output.
#' @param type Character. The default, "joint" creates a grid of plots, one plot per time value, showing training data incremented with the predictions. ""SD" creates a grid of plots showing standard deviations for predictions. "trainingData" and "predictions" produce a similar grid with only the training data and predictions, respectively. "pars" produces plots of the marginal posterior densities for the parameters and hyperparameters.
#' @param polygonsToOverlay Optional SpatialPolygons object. Those polygons, e.g. regional boundaries, are added to each prediction map.
#' @param ... Arguments for the graphics engine; 'width' and 'height' should be especially helpful.
#' @param control List of control options, cf. ?plot.control.
#'
#' @details The function produces a grid of rasters by default, one per time value. If your data are not spaced regularly, you can plot the result as points instead by setting `control = plot.control(plotRaster = FALSE)`. If you have data collected over more than a few days, selecting the days to plot with the `timeValues` option in control is recommended.
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
plot.INLAMRA <- function(x, filename = NULL, type = c("joint", "training", "predictions", "SD", "marginals"), polygonsToOverlay = NULL, control = NULL, ...) {
  if (is.null(control)) {
    control <- plot.control()
  } else {
    control <- do.call(plot.control, args = control)
  }
  type <- type[[1]]

  .checkPlotInputConsistency(x, type)

  if (control$trim == 1) {
    x$predMoments$Mean <- replace(x$predMoments$Mean, which(x$predMoments$Mean > max(x$data$spObject@data$y)), max(x$data$spObject@data$y))
    x$predMoments$Mean <- replace(x$predMoments$Mean, which(x$predMoments$Mean < min(x$data$spObject@data$y)), min(x$data$spObject@data$y))
  } else if (control$trim == 2) {
    for (timeValue in unique(x$predData$time)) {
      if (timeValue %in% x$data$time) {
        minValueThatDay <- min(x$data$spObject@data$y[x$data$time == timeValue])
        maxValueThatDay <- max(x$data$spObject@data$y[x$data$time == timeValue])
        x$predMoments$Mean <- replace(x$predMoments$Mean, which((x$predMoments$Mean > maxValueThatDay) & (x$predData$time == timeValue)), maxValueThatDay)
        x$predMoments$Mean <- replace(x$predMoments$Mean, which((x$predMoments$Mean < minValueThatDay) & (x$predData$time == timeValue)), minValueThatDay)
      }
    }
  }

  if (control$plotRaster & (type != "marginals")) {
    plottedObjects <- .plotRasters(x, control, filename = filename, type = type, polygonsToOverlay = polygonsToOverlay, ...)
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

.plotRasters <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type, ...) {

  rasterListPerTimeUnit <- .rasterizeTrainingTestJointSD(INLAMRAoutput = INLAMRAoutput, control = control)

  rastersToPlot <- lapply(rasterListPerTimeUnit, function(x) x[[type]])
  names(rastersToPlot) <- names(rasterListPerTimeUnit)

  numNonEmptyPlots <- sum(sapply(rastersToPlot, FUN = function(x) any(!is.na(raster::values(x)))))
  numPlotsPerLine <- ceiling(sqrt(numNonEmptyPlots))
  if (!is.null(filename)) {
    control$graphicsEngine(filename = filename, ...)
  }
  layout(matrix(1:(numPlotsPerLine^2), nrow = numPlotsPerLine, ncol = numPlotsPerLine, byrow = TRUE))

  valuesRanges <- sapply(rastersToPlot, function(x) range(raster::values(x), na.rm = TRUE))
  boundaries <- c(min(valuesRanges[1, ]) - 1e-300, max(valuesRanges[2, ]) + 1e-300)
  colourBreaks <- seq(from = boundaries[[1]], to = boundaries[[2]],                      length.out = ifelse(is.null(control$controlForRasterColourScale$breaks), 10, control$controlForRasterColourScale$breaks))
  ecol <- do.call(mapmisc::colourScale, args = c(list(x = raster::values(rastersToPlot[[1]])), control$controlForRasterColourScale))
  for (i in seq_along(rastersToPlot)) {
    if (any(!is.na(raster::values(rastersToPlot[[i]])))) {
      if (!("main" %in% names(control$controlForRasterPlot))) {
        do.call(raster::plot, args = c(list(x = rastersToPlot[[i]], col = ecol$col, breaks = colourBreaks, legend = FALSE, main = names(rastersToPlot)[[i]]), control$controlForRasterPlot)) # I'm not using the breaks set by colourScale because I want the breaks to be the exact same across all graphs, which are supposed to express similar quantities on the same scale.
      } else {
        do.call(raster::plot, args = c(list(x = rastersToPlot[[i]], col = ecol$col, breaks = colourBreaks, legend = FALSE), control$controlForRasterPlot))
      }
      # do.call(raster::plot, args = c(list(x = rastersToPlot[[i]], legend.only = TRUE, zlim = colorRange), control$controlForRasterLegend))

      if (!is.null(polygonsToOverlay)) {
        subPolygons <- raster::intersect(polygonsToOverlay, raster::extent(rastersToPlot[[i]]))
        raster::plot(subPolygons, add = TRUE)
      }

      if (!is.na(raster::crs(rastersToPlot[[i]]))) {
        do.call(mapmisc::scaleBar, args = c(list(crs = raster::crs(rastersToPlot[[i]])), control$controlForScaleBar))
      }

      do.call(mapmisc::legendBreaks, args = c(list(breaks = ecol), control$controlForRasterLegend))
    }
  }
  if (!is.null(filename)) {
    dev.off()
  }
  rastersToPlot
}

# Only applied to time points where predictions were required.

.rasterizeTrainingTestJointSD <- function(INLAMRAoutput, control) {
  combinedCoordMat <- rbind(INLAMRAoutput$data$spObject@coords, INLAMRAoutput$predData$spObject@coords)
  mapExtentAndRasterDims <- .getExtentAndRasterSizes(coordMat = combinedCoordMat, resolutionInMeters = control$resolutionInMeters)
  landRaster <- raster::raster(nrows = mapExtentAndRasterDims$rasterDims$nRows, ncols = mapExtentAndRasterDims$rasterDims$nCols, ext = mapExtentAndRasterDims$extentObj, crs = sp::proj4string(INLAMRAoutput$data$spObject))

  timeValues <- sort(unique(c(INLAMRAoutput$data$time, INLAMRAoutput$predData$time)))
  if (!is.null(control$timesToPlot)) {
    if (!any(control$timesToPlot %in% timeValues)) {
      stop("The time points inputted for plotting through timesToPlot are not found in the data. Make sure they were inputted correctly. \n")
    }
    timeValues <- sort(control$timesToPlot)
  }
  if (length(timeValues) > 14) {
    warning("You will be getting rasters for more than 14 time points. Result might look bad. Either specify a subset of time points with control$timesToPlot or, if you have spatiotemporal data collected at irregular intervals, set control$plotRaster to FALSE.")
  }
  funForRasterList <- function(timePoint) {
    trainingDataIndices <- INLAMRAoutput$data$time == timePoint
    testDataIndices <- INLAMRAoutput$predData$time == timePoint

    rasterList <- list()
    rasterList$SD <- rasterList$predictions <- rasterList$training <- landRaster

    if (any(testDataIndices)) {
      rasterList$SD <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$SD[testDataIndices])
      rasterList$predictions <- raster::rasterize(x = INLAMRAoutput$predData$spObject@coords[testDataIndices, ], y = landRaster, field = INLAMRAoutput$predMoments$Mean[testDataIndices])
    }
    if (any(trainingDataIndices)) {
      rasterList$training <- rasterList$joint <- raster::rasterize(x = INLAMRAoutput$data$spObject@coords[trainingDataIndices, ], y = landRaster, field = INLAMRAoutput$data$spObject@data[trainingDataIndices , "y"])
      if (any(testDataIndices)) {
        jointCoordinates <- unname(rbind(INLAMRAoutput$predData$spObject@coords[testDataIndices, ], INLAMRAoutput$data$spObject@coords[trainingDataIndices, ]))
        dataObject <- data.frame(y = unname(c(INLAMRAoutput$predMoments$Mean[testDataIndices], INLAMRAoutput$data$spObject@data[trainingDataIndices, "y"])))
        pointsDataFrame <- sp::SpatialPointsDataFrame(coords = jointCoordinates, data = dataObject)
        rasterList$joint <- raster::rasterize(x = pointsDataFrame, y = landRaster, field = "y")
      }
    } else {
      rasterList$joint <- rasterList$predictions
    }
    rasterList
  }
  funForRasterListOutputs <- lapply(timeValues, funForRasterList)
  namesForElements <- paste("Time =", as.character(timeValues))
  names(funForRasterListOutputs) <- namesForElements
  funForRasterListOutputs
}

.getExtentAndRasterSizes <- function(coordMat, resolutionInMeters) {
  rangeByColumn <- apply(coordMat, MARGIN = 2, range)
  cornerCoordinates <- rbind(
    rangeByColumn[1, ],
    replace(rangeByColumn[1, ], 1, rangeByColumn[2, 1]),
    rangeByColumn[2, ],
    replace(rangeByColumn[1, ], 2, rangeByColumn[2, 2])
  )
  lonDistInMeters <- geosphere::distHaversine(p1 = cornerCoordinates[1, ], p2 = cornerCoordinates[2, ])
  latDistInMeters <- geosphere::distHaversine(p1 = cornerCoordinates[1, ], p2 = cornerCoordinates[4, ])
  # The + 1 is to account for the fact that the corners of the zone identified fall in the middle of raster tiles.
  rasterSizes <- list(
    nCols = round(lonDistInMeters/resolutionInMeters) + 1,
    nRows = round(latDistInMeters/resolutionInMeters) + 1
  )
  lonPadding <- (cornerCoordinates[2, 1] - cornerCoordinates[1, 1])/(2*(rasterSizes$nCols - 1))
  latPadding <- (cornerCoordinates[4, 2] - cornerCoordinates[1, 2])/(2*(rasterSizes$nRows - 1))
  extentBox <- rangeByColumn
  extentBox[1, ] <- extentBox[1, ] - c(lonPadding, latPadding)
  extentBox[2, ] <- extentBox[2, ] + c(lonPadding, latPadding)

  list(extentObj = raster::extent(t(extentBox)), rasterDims = rasterSizes)
}

.plotPoints <- function(INLAMRAoutput, control, filename, polygonsToOverlay, type, ...) {
  pointsPerDayList <- .getPointsPerDay(INLAMRAoutput = INLAMRAoutput, control = control)
  pointsToPlotList <- lapply(pointsPerDayList, function(x) x[[type]])

  pointsForRange <- lapply(pointsToPlotList, function(x) x[["joint"]])
  pointColourRanges <- sapply(pointsForRange, function(x) range(x@data$y))
  colorRange <- c(min(pointColourRanges[1, ]), max(pointColourRanges[2, ]))

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

#' Control parameters for plot.INLAMRA
#'
#' Control parameters for the `control` argument of plot.
#'
#' @param trim Either 0, 1, or 2. Determines if predicted values should be trimmed before plotting to fit within the observation range. 0 (default): no trimming, 1: Values are trimmed to fall within the range of all observations, 2: Values are trimmed to fall within the range of observations on the corresponding day (if available).
#' @param plotRaster Logical. Should rasters (TRUE) or points (FALSE) be plotted? Ignored if type is "marginals". Points may be preferrable if data were not collected on a regular spatial grid.
#' @param graphicsEngine Function, e.g. jpeg, png, postscript. Corresponds to the graphics engine used to produce the plots,
#' @param resolutionInMeters Numeric value indicating the width, in meters, of each (square) cell in the raster. Ignored if plotRaster is FALSE.
#' @param timesToPlot Vector of dates in the same format as that used in the data used to fit INLAMRA. Only maps corresponding to those days will be plotted.
#' @param controlForScaleBar List of control parameters for the scale bars in the graphs, cf. ?mapmisc::scaleBar.
#' @param controlForRasterPlot List of control parameters for the raster plotting function, cf. ?raster::plot.
#' @param controlForRasterLegend List of control parameters for the raster legend, cf. ?mapmisc::legendBreaks.
#' @param controlForRasterColourScale List of control parameters for the raster colour scale, cf. ?mapmisc::colourScale.
#'
#' @details The function need not be called explicitly: it's just a convenient way to see/set control parameters.
#'
#' @return A list of control parameters.
#'
#' @examples
#' \dontrun{
#' }
#' @export
#'

plot.control <- function(trim = 0, plotRaster = TRUE, graphicsEngine = jpeg, resolutionInMeters = NULL, timesToPlot = NULL, controlForScaleBar = NULL, controlForRasterPlot = NULL, controlForRasterLegend = NULL, controlForRasterColourScale = NULL) {
  list(trim = trim, plotRaster = plotRaster, graphicsEngine = graphicsEngine, controlForScaleBar = controlForScaleBar, controlForRasterPlot = controlForRasterPlot, controlForRasterLegend = controlForRasterLegend, controlForRasterColourScale = controlForRasterColourScale, resolutionInMeters = resolutionInMeters, timesToPlot = timesToPlot)
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
