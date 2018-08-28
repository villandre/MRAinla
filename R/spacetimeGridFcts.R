#' Constructs an irregular grid on a specified spatiotemporal range.
#'
#' Constructing a spatiotemporal grid will be necessary prior to fitting the spatiotemporal MRA
#' proposed by Villandre et al.
#'
#' @param lonBreaks vector indicating where the longitude boudaries are located
#' @param latBreaks vector indicating where the latitude boudaries are located
#' @param timeBreaks vector indicating where the temporal boudaries are located
#' @param numKnots either the number of knots in all elements of the grids, or a vector of size equal to
#' the length of gridRasterList giving the number of knots for each resolution, or a list of vectors giving
#' the number of knots in each region of the grid for each resolution
#' @param hyperPriorFunList named list of functions with one argument specifying the hyperprior distributions
#'
#' @details A 3D irregular grid can be represented parsemoniously by a list of 3 vectors, one #' for each dimension. Under this scheme, spatiotemporal regions at a given resolution are
#' like bricks. This is a constructor function.
#'
#' @return A Spacetimegrid object.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

gridConstructor <- function(lonBreaks, latBreaks, timeBreaks, lonExtent, latExtent, timeExtent) {
  grids <- mapply(breaks = list(lonBreaks, latBreaks, timeBreaks), extent = list(lonExtent, latExtent, timeExtent), FUN = function(breaks, extent) {
    list(breaks = breaks, extent = extent)
  }, SIMPLIFY = FALSE)
  names(grids) <- c("longitude", "latitude", "time")
  grids$midpoints <- .deriveMidpoints(grids)
  class(grids) <- "Spacetimegrid"
  grids
}

#' A custom print function for the Spacetimegrid class.
#'
#' Prints a Spacetimegrid object.
#'
#' @param x Spacetimegrid object
#'
#' @details
#'
#' @return Print function, doesn't return anything.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

print.Spacetimegrid <- function(x) {
  cat("Spacetime grid object \n")
  cat("Depth: ", gridDepth <- .getM(x), "\n")
}

#' A custom plot function for the Spacetimegrid class.
#'
#' Plots the spatial grid within the Spacetimegrid object, plus, optionally, observations and knots.
#'
#' @param x Spacetimegrid object
#' @param observationsAsPoints SpatialPoints or SpatialPointsDataFrame object containing the
#' observations
#' @param knotsAsPoints SpatialPoints or SpatialPointsDataFrame object containing the
#' knots
#'
#' @details The function plots a 2D projection of the grid, flattening the dataset with respect to time.
#'
#' @return Plot function, doesn't return anything worthwhile.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

plot.Spacetimegrid <- function(x, observationsAsPoints, knotsAsPoints) {
  combinedSP <- NULL
  plotColours <- NULL
  if (!is.null(observationsAsPoints)) {
    observSP <- SpatialPointsDataFrame(observationsAsPoints, data = data.frame(identity = rep("red", nrow(observationsAsPoints@coords))))
    combinedSP <- observSP
    plotColours <- combinedSP@data$identity
  }
  if (!is.null(knotsAsPoints)) {
    knotSP <- SpatialPointsDataFrame(knotsAsPoints, data = data.frame(identity = rep("green", nrow(knotsAsPoints@coords))))
    combinedSP <- raster::bind(combinedSP, knotSP)
    plotColours <- combinedSP@data$identity
  }
  plot(combinedSP, col = as.character(plotColours), xlab = "Longitude", ylab = "Latitude", xlim = x$longitude$extent, ylim = x$latitude$extent, cex = 0.8, pch = 18)
  abline(v = c(x$longitude$extent[1], x$longitude$breaks, x$longitude$extent[2]), lty = 3, col = "blue4")
  abline(h = c(x$latitude$extent[1], x$latitude$breaks, x$latitude$extent[2]), lty = 3, col = "blue4")
}

# observations are coded as SpatialPointsDataFrame with 3D coordinates, longitude-latitude-time
# In each dimension, we have ranges defined [.,.), i.e. closed on the left, open on the right.

.getSpacetimeDim <- function(spacetimeObj, dimension = c("longitude", "latitude", "time")) {
  dimension <- dimension[[1]]
  if (dimension == "longitude") {
    return(spacetimeObj@sp@coords[ , 1])
  }
  if (dimension == "latitude") {
    return(spacetimeObj@sp@coords[ , 2])
  }
  return(index(spacetimeObj))
}

#' Add breaks in space-time grid.
#'
#' This function is used to refine a resolution grid.
#'
#' @param spacetimegridObj Spacetimegrid object
#' @param dimension the name of the dimension where breaks will be added, either "longitude",
#' "latitude", or "time"
#' @param breaks vector indicating where the new breakpoints should be located.
#'
#' @details This function is meant to create an embedded resolution.
#'
#' @return A Spacetimegrid object with the new breaks.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export
addBreaks <- function(spacetimegridObj, dimension = c("longitude", "latitude", "time"), breaks) {
  dimension <- dimension[[1]]
  spacetimegridObj[[dimension]]$breaks <- sort(c(spacetimegridObj[[dimension]]$breaks, breaks))
  newMidpoints <- .deriveMidpoints(spacetimegridObj)
  spacetimegridObj$midpoints <- newMidpoints
  spacetimegridObj
}

.sameGridSection <- function(x, y, spacetimegridObj) {
  sameSection <- FALSE
  testResultTime <- identical(match(TRUE, spacetimegridObj[["time"]]$breaks > index(x)), match(TRUE, spacetimegridObj[["time"]]$breaks > index(y)))
  if (testResultTime) {
    testResultsSpace <- sapply(c("longitude", "latitude"), FUN = function(dimensionName) {
      xUpperPos <- match(TRUE, spacetimegridObj[[dimensionName]]$breaks > .getSpacetimeDim(x, dimensionName))
      yUpperPos <- match(TRUE, spacetimegridObj[[dimensionName]]$breaks > .getSpacetimeDim(y, dimensionName))
      identical(xUpperPos, yUpperPos)
    })
    if (all(testResultsSpace)) {
      sameSection <- TRUE
    }
  }
  sameSection
}

.deriveMidpoints <- function(partialGridObj) {
  midpointsPerDim <- lapply(partialGridObj[c("longitude", "latitude", "time")], FUN = function(dimension) {
    extendedBreaks <- c(dimension$extent[[1]], dimension$breaks, dimension$extent[[2]])
    if ("POSIXct" %in% class(extendedBreaks)) {
      return(as.POSIXct((as.numeric(head(extendedBreaks, n = -1)) + as.numeric(tail(extendedBreaks, n = -1)))/2, origin = "1970-01-01"))
    }
    (head(extendedBreaks, n = -1) + tail(extendedBreaks, n = -1))/2
  })
  allCombinations <- expand.grid(seq_along(midpointsPerDim[[1]]), seq_along(midpointsPerDim[[2]]), seq_along(midpointsPerDim[[3]]))
  allPointsSpaceAndTime <- lapply(1:nrow(allCombinations), FUN = function(coordTrioIndex) {
    aPoint <- sp::SpatialPoints(matrix(c(midpointsPerDim[[1]][allCombinations[coordTrioIndex, 1]], midpointsPerDim[[2]][allCombinations[coordTrioIndex, 2]]), nrow = 1, ncol = 2))
    timeCoord <- midpointsPerDim[[3]][allCombinations[coordTrioIndex, 3]]
  list(space = aPoint, time = timeCoord)
  })
  pointsList <- lapply(allPointsSpaceAndTime, FUN = function(x) x$space)
  if (length(pointsList) > 1) {
    pointsMerged <- do.call(raster::bind, pointsList)
  } else {
    pointsMerged <- pointsList[[1]]
  }
  STI(sp = pointsMerged, time = do.call("c", lapply(allPointsSpaceAndTime, function(x) x$time)))
}

# gridConstructorAlt <- function(spacetimegridObj = NULL, lonBreaks, latBreaks, timeBreaks, lonExtent, latExtent, timeExtent) {
#   if (is.null(spacetimeGridObj)) {
#     extendedBreakList <- list(longitude = c(min(lonBreaks), lonBreaks, max(lonBreaks)), latitude = c(min(latBreaks), latBreaks, max(latBreaks)), time = c(min(timeBreaks), timeBreaks, max(timeBreaks)))
#     allRanges <- lapply(extendedBreakList, FUN = function(breaks) {
#       cbind(head(extendedLon, n = -1), tail(extendedLon, n = -1))
#     })
#     names(allRanges) <- names(extendedBreakList)
#     numIntervals <- sapply(allRanges, nrow)
#     combinations <- expand.grid(1:numIntervals[[1]], 1:numIntervals[[2]], 1:numIntervals[[3]])
#     colnames(combinations) <- names(allRanges)
#     allBricks <- lapply(1:nrow(combinations), FUN = function(combIndex) {
#       combination <- combinations[combIndex,, drop = FALSE]
#       gridlayerObject <- list(lonExtent = allRanges$longitude[combIndex$longitude], latExtent = allRanges$latitude[combIndex$latitude], timeExtent = allRanges$time[combIndex$time])
#       class(gridlayerObject) <- "SpacetimegridLayer"
#       gridlayerObject
#     })
#     return(allBricks)
#   }
#   updatedGrid <- gridRecursion(spacetimegridObj, lonBreaks, latBreaks, timeBreaks)
# }

# gridRecursion <- function(spacetimegridObj, lonBreaks, latBreaks, timeBreaks) {
#   newBreaksList <- list(longitude = lonBreaks, latitude = latBreaks, time = timeBreaks)
#   if ("SpacetimegridLayer" %in% class(spacetimegridObj[[1]])) {
#     gridRecursion(spacetimegridObj, lonBreaks, latBreaks, timeBreaks)
#   }
#   ### ADD RECURSIVE COMPONENT ############
#   lapply(spacetimegridObj, function(spacetimeBrick) {
#     newBreaksForBrick <- lapply(names(newBreaksList), FUN = function(dimName) {
#       test1 <- newBreaksList[[dimName]] < max(get(spacetimeBrick, dimName))
#       test2 <- newBreaksList[[dimName]] > min(get(spacetimeBrick, dimName))
#       newBreaks <- NULL
#       if (any(test1 & test2)) {
#         breakToAddIndices <- which(test1 & test2)
#         newBreaks <- newBreaksList[breakToAddIndices]
#       }
#       newBreaks
#     })
#     names(newBreaksForBrick) <- names(newBreaksList)
#     if (all(is.null(sapply(newBreaksForBrick, function(x) x)))) {
#       return(spacetimeBrick)
#     }
#     newBricks <- splitBrick(lonBreaks = newBreaksForBrick$longitude, latBreaks = newBreaksForBrick$latitude, timeBreaks = newBreaksForBrick$time)
#     return(newBricks)
#   })
# }

spacetimebrickConstructor <- function(lonExtent, latExtent, timeExtent, parentBrick = NULL, observations = NULL) {
  if (is.null(parentBrick)) {
    parentBrick <- .GlobalEnv
  }
  if (!is.null(observations)) {
    observations <- subset(observations, lonExtent = lonExtent, latExtent = latExtent, timeExtent = timeExtent)
  }
  brickEnvironment <- new.env(parent = parentBrick)
  assign(x = "dimensions", value = list(longitude = lonExtent, latitude = latExtent, time = timeExtent), envir = brickEnvironment)
  assign(x = "knotPositions", value =  NULL, envir = brickEnvironment)
  assign(x = "observations", value = observations, envir = brickEnvironment)
  assign(x = "childBricks", value = NULL, envir = brickEnvironment)
  if (!identical(parentBrick, .GlobalEnv)) {
    childAddresses <- get(x = "childBricks", envir = parentBrick)
    incAddresses <- c(childAddresses, brickEnvironment)
    assign(x = "childBricks", value = incAddresses, envir = parentBrick)
  }
  class(brickEnvironment) <- "Spacetimebrick"
  brickEnvironment
}

# Function is used to define an initial grid, nested by default in the resolution 0 image. We can use the same function to nest a grid in any spacetimebrick.
# Initial breaks should include the lower and upper limit of each dimension.
# When a parent brick is designated, it is split based on the breaks in lonBreaks, latBreaks and timeBreaks. In that situation, the function keeps only breaks defined in the space-time region of the brick.

SpacetimegridConstructor <- function(parentBrick = NULL, lonBreaks, latBreaks, timeBreaks, observations = NULL) {
  if (is.null(parentBrick)) { # Create a trivial brick which all generated bricks will refer to.
    parentBrick <- spacetimebrickConstructor(lonExtent = range(lonBreaks), latExtent = range(latBreaks), timeExtent = range(timeBreaks), observations = observations)
  }
  breakList <- list(longitude = lonBreaks, latitude = latBreaks, time = timeBreaks)
  correctedBreakList <- lapply(c("longitude", "latitude", "time"), FUN = function(dimName) {
    extent <- get(x = "dimensions", envir = parentBrick)[[dimName]]
    proposedBreaks <- breakList[[dimName]]

    keepBreakIndices <- (proposedBreaks < max(extent)) & (proposedBreaks > min(extent))
    remainingBreaks <- c(min(extent), proposedBreaks[keepBreakIndices], max(extent))
  })
  allRanges <- lapply(correctedBreakList, FUN = function(breaks) {
    cbind(head(breaks, n = -1), tail(breaks, n = -1))
  })
  names(allRanges) <- names(breakList)
  numIntervals <- sapply(allRanges, nrow)
  combinations <- expand.grid(1:numIntervals[[1]], 1:numIntervals[[2]], 1:numIntervals[[3]])
  colnames(combinations) <- names(allRanges)
  lapply(1:nrow(combinations), FUN = function(combIndex) {
    combination <- combinations[combIndex,, drop = FALSE]
    spacetimebrickConstructor(lonExtent = allRanges$longitude[combination$longitude,], latExtent = allRanges$latitude[combination$latitude,], timeExtent = as.POSIXct(allRanges$time[combination$time,], origin = '1970-01-01'), parentBrick = parentBrick, observations = parentBrick$observations)
  })
  topAddress <- getTopEnvirAddress(parentBrick)
  class(topAddress) <- "Spacetimegrid"
  topAddress
}

getTopEnvirAddress <- function(nestedEnvir) {
  parentEnvir <- parent.env(nestedEnvir)
  if (identical(parentEnvir, .GlobalEnv)) {
    return(nestedEnvir)
  }
  getTopEnvirAddress(parentEnvir)
}

.getM <- function(spacetimegridObj) {
  counter <- 0
  currentAddress <- spacetimegridObj
  while (!is.null(currentAddress$childBricks)) {
    counter <- counter + 1
    currentAddress <- currentAddress$childBricks[[1]]
  }
  counter
}

subset.STI <- function(x, latExtent, lonExtent, timeExtent) {
  timeIndices <- (index(x) <= max(timeExtent)) & (index(x) > min(timeExtent))
  lonIndices <- (x@sp@coords[, 1] <= max(lonExtent)) & (x@sp@coords[, 1] > min(lonExtent))
  latIndices <- (x@sp@coords[, 2] <= max(latExtent)) & (x@sp@coords[, 2] > min(latExtent))
  x[which(timeIndices & lonIndices & latIndices)]
}

addLayer <- function(spacetimegridObj, latBreaks, lonBreaks, timeBreaks) {
  tipAddresses <- .tipAddresses(spacetimegridObj)
  lapply(tipAddresses, FUN = function(tipAddress) {
    SpacetimegridConstructor(parentBrick = tipAddress, lonBreaks = lonBreaks, latBreaks = latBreaks, timeBreaks = timeBreaks)
    NULL
  })
  cat("Done. Beware of side-effects! \n")
}

.tipAddresses <- function(spacetimebrickObj) {
  if (!is.null(spacetimebrickObj$childBricks)) {
    return(lapply(spacetimebrickObj$childBricks, FUN = .tipAddresses))
  }
  spacetimebrickObj
}
