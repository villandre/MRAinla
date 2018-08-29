#' Constructs an irregular grid on a specified spatiotemporal range.
#'
#' Constructing a spatiotemporal grid will be necessary prior to fitting the spatiotemporal MRA
#' proposed by Villandre et al.
#'
#' @param lonBreaks vector indicating where the longitude boudaries and breaks are located
#' @param latBreaks vector indicating where the latitude boudaries and breaks are located
#' @param timeBreaks vector indicating where the temporal boudaries and breaks are located
#' @param observations optional STI object specifying where the observations are located.
#' @param parentBrick (used for nesting grids, don't touch this) optional spacetimebrick object.
#'
#' @details # This function is used to define an initial grid, nested by default in the
#' resolution 0 image. Initial breaks should include the lower and upper limit of each dimension. It is also possible to define the grid within an existing space-time brick
#' object, called the 'parent brick'.#'
#' @return A Spacetimegrid object.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

SpacetimegridConstructor <- function(lonBreaks, latBreaks, timeBreaks, observations = NULL, parentBrick = NULL) {
  breakList <- list(longitude = lonBreaks, latitude = latBreaks, time = timeBreaks)
  if (is.null(parentBrick)) { # Create a trivial brick which all generated bricks will refer to.
    parentBrick <- .spacetimebrickConstructor(lonExtent = range(lonBreaks), latExtent = range(latBreaks), timeExtent = range(timeBreaks), observations = observations)
    parentBrick$breaks <- list(breakList)
  }
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
    .spacetimebrickConstructor(lonExtent = allRanges$longitude[combination$longitude,], latExtent = allRanges$latitude[combination$latitude,], timeExtent = as.POSIXct(allRanges$time[combination$time,], origin = '1970-01-01'), parentBrick = parentBrick, observations = parentBrick$observations)
  })
  topAddress <- .getTopEnvirAddress(parentBrick)
  class(topAddress) <- "Spacetimegrid"
  topAddress
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
 # TO DO
}

#' Add breaks in space-time grid.
#'
#' This function is used to refine a resolution grid.
#'
#' @param spacetimegridObj Spacetimegrid object
#' @param latBreaks vector indicating where the new latitude breakpoints should be located.
#' @param lonBreaks vector indicating where the new longitude breakpoints should be located.
#' @param timeBreaks vector indicating where the new time breakpoints should be located.
#'
#' @details This function creates an embedded resolution.
#'
#' @return A Spacetimegrid object with a new layer.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

addLayer <- function(spacetimegridObj, latBreaks, lonBreaks, timeBreaks) {
  newBreaks <- list(longitude = sort(c(tail(spacetimegridObj$breaks, n = 1)[[1]]$longitude, lonBreaks)), latitude = sort(c(tail(spacetimegridObj$breaks, n = 1)[[1]]$latitude, latBreaks)), time = sort(c(tail(spacetimegridObj$breaks, n = 1)[[1]]$time, timeBreaks)))
  spacetimegridObj$breaks <- c(spacetimegridObj$breaks, list(newBreaks))
  tipAddresses <- .tipAddresses(spacetimegridObj)
  lapply(tipAddresses, FUN = function(tipAddress) {
    SpacetimegridConstructor(parentBrick = tipAddress, lonBreaks = lonBreaks, latBreaks = latBreaks, timeBreaks = timeBreaks)
    NULL
  })
  cat("Done. Beware of side-effects! \n")
}

#' Add breaks in space-time grid.
#'
#' This function is used to refine a resolution grid.
#'
#' @param spacetimegridObj Spacetimegrid object
#' @param latBreaks vector indicating where the new latitude breakpoints should be located.
#' @param lonBreaks vector indicating where the new longitude breakpoints should be located.
#' @param timeBreaks vector indicating where the new time breakpoints should be located.
#'
#' @details This function creates an embedded resolution.
#'
#' @return A Spacetimegrid object with a new layer.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

getLayer <- function(spacetimeGridObj, m = 0) {
  if (m == 0) {
    return(spacetimeGridObj)
  }
  nestedResults <- .dive(x = spacetimeGridObj, m = m-1)
  unlist(nestedResults, recursive = TRUE)
}

.sameGridSection <- function(x, y, m, spacetimegridObj) {
  if (m == 0) {
    return(TRUE)
  }
  testResultTime <- identical(match(TRUE, spacetimegridObj$breaks[[m]][["time"]] > .getSpacetimeDim(x, "time")), match(TRUE, spacetimegridObj$breaks[[m]][["time"]] > .getSpacetimeDim(y, "time")))
  if (!testResultTime) {
    return(FALSE)
  }
  testResultsSpace <- sapply(c("longitude", "latitude"), FUN = function(dimensionName) {
    xUpperPos <- match(TRUE, spacetimegridObj$breaks[[m]][[dimensionName]] > .getSpacetimeDim(x, dimensionName))
    yUpperPos <- match(TRUE, spacetimegridObj$breaks[[m]][[dimensionName]] > .getSpacetimeDim(y, dimensionName))
    identical(xUpperPos, yUpperPos)
  })
  all(testResultsSpace)
}

.spacetimebrickConstructor <- function(lonExtent, latExtent, timeExtent, parentBrick = NULL, observations = NULL) {
  if (is.null(parentBrick)) {
    parentBrick <- emptyenv()
  }
  if (!is.null(observations)) {
    observations <- subset(observations, lonExtent = lonExtent, latExtent = latExtent, timeExtent = timeExtent)
  }
  brickEnvironment <- new.env(parent = parentBrick)
  assign(x = "dimensions", value = list(longitude = lonExtent, latitude = latExtent, time = timeExtent), envir = brickEnvironment)
  assign(x = "knotPositions", value =  NULL, envir = brickEnvironment)
  assign(x = "observations", value = observations, envir = brickEnvironment)
  assign(x = "childBricks", value = NULL, envir = brickEnvironment)
  if (!identical(parentBrick, emptyenv())) {
    childAddresses <- get(x = "childBricks", envir = parentBrick)
    incAddresses <- c(childAddresses, brickEnvironment)
    assign(x = "childBricks", value = incAddresses, envir = parentBrick)
  }
  class(brickEnvironment) <- "Spacetimebrick"
  brickEnvironment
}



.getTopEnvirAddress <- function(nestedEnvir) {
  parentEnvir <- parent.env(nestedEnvir)
  if (identical(parentEnvir, emptyenv())) {
    return(nestedEnvir)
  }
  .getTopEnvirAddress(parentEnvir)
}

.getM <- function(spacetimegridObj) {
  length(spacetimegridObj$breaks)
}

subset.STI <- function(x, latExtent, lonExtent, timeExtent) {
  timeIndices <- (index(x) <= max(timeExtent)) & (index(x) > min(timeExtent))
  lonIndices <- (x@sp@coords[, 1] <= max(lonExtent)) & (x@sp@coords[, 1] > min(lonExtent))
  latIndices <- (x@sp@coords[, 2] <= max(latExtent)) & (x@sp@coords[, 2] > min(latExtent))
  x[which(timeIndices & lonIndices & latIndices)]
}

.tipAddresses <- function(spacetimebrickObj) {
  if (!is.null(spacetimebrickObj$childBricks)) {
    return(lapply(spacetimebrickObj$childBricks, FUN = .tipAddresses))
  }
  spacetimebrickObj
}

.dive <- function(x, m) {
  if (m == 0) {
    return(x$childBricks)
  }
  lapply(x$childBricks, .dive, m = m-1)
}

#' Specify knots at each resolution.
#'
#' This function lets a user specify the knot positions at every resolution. If desired, it will automatically place an arbitrary number of equidistant knots in each section of the grid.
#'
#' @param spacetimegridObj Spacetimegrid object
#' @param knotsList list of knot positions in spacetime format. If it has length 1, it will be assumed to be the same for all dimensions
#' @param r scalar or vector indicating the number of knots in each section of the grid at each resolution. Ignored if knotsList is specified. If it is a scalar, it is assumed that the number of knots is the same for each section of the grid at every resolution.
#'
#' @details This function operates with side-effects! It modifies spacetimegridObj, specifying the knotPositions component of each Spacetimebrick object.
#'
#' @return Quietly, an object of class Spacetimegrid. Note that it directly modifies spacetimegridObj. In other words, it relies on side-effects.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

addKnots <- function(spacetimegridObj, knotsList = NULL, r = 10) {
  if (("STI" %in% class(knotsList))) { # If knotsList is NULL, class(knotsList) will return "NULL" (the class of NULL if "NULL", not inexistent.)
    knotsList <- replicate(.getM(spacetimegridObj)+1, expr = knotsList, simplify = FALSE) # In this situation, the knots are the same at every resolution. Might not be the best idea...
  }
  if (is.list(knotsList) & (length(knotsList) == 1)) {
    knotsList <- replicate(.getM(spacetimegridObj)+1, expr = knotsList[[1]], simplify = FALSE)
  }
  .populateKnots(spacetimegridObj, knotsList, r = r)

  spacetimegridObj
}

.populateKnots <- function(spacetimegridObj, knotsList, r, tuningPara = 0.1) {
  if (is.null(knotsList)) {
    currentR <- r[[1]]
    lowerTime <- min(spacetimegridObj$dimensions$time)
    upperTime <- floor(max(spacetimegridObj$dimensions$time)-1)
    timeCoords <- rep(c(lowerTime, upperTime), length.out = currentR)
    numLines <- ceiling(r/2)
    baseLonCoords <- matrix(c(1-tuningPara, tuningPara, tuningPara, 1-tuningPara,2),2)%*%spacetimegridObj$dimensions$longitude
    lonCoords <- rep(baseLonCoords, length.out = currentR)
    latStepSize <- diff(spacetimegridObj$dimensions$latitude)
    baseLatCoords <- seq(from = min(spacetimegridObj$dimensions$latitude)+latStepSize, to = max(spacetimegridObj$dimensions$latitude)-latStepSize, each = numLines, length.out = currentR)
    knotsList <- list(STI(sp = SpatialPoints(cbind(lonCoords, latCoords), time = timeCoords)))
  }
  spacetimegridObj$knotPositions <- knotsList[[1]]
  if (!is.null(spacetimegridObj$childBricks)) {
    lapply(spacetimegridObj$childBricks, .populateKnots, knotsList = tail(knotsList, n = -1), r = tail(r, n = -1), tuningPara = tuningPara)
  }
  invisible()
}
