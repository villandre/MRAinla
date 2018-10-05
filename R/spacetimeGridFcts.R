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

.SpacetimegridConstructor <- function(lonBreaks, latBreaks, timeBreaks, observations = NULL, parentBrick = NULL) {
  breakList <- list(longitude = lonBreaks, latitude = latBreaks, time = timeBreaks)

  gridCreationFlag <- FALSE

  if (is.null(parentBrick)) { # Create a trivial brick which all generated bricks will refer to.
    gridCreationFlag <- TRUE
    parentBrick <- .SpacetimebrickConstructor(lonExtent = range(lonBreaks), latExtent = range(latBreaks), timeExtent = range(timeBreaks), observations = observations)
    parentBrick$breaks <- list(breakList)
  }

  correctedBreakList <- lapply(c("longitude", "latitude", "time"), FUN = function(dimName) {
    extent <- parentBrick$dimensions[[dimName]]
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
    combination <- combinations[combIndex, , drop = FALSE]
    .SpacetimebrickConstructor(lonExtent = allRanges$longitude[combination$longitude,], latExtent = allRanges$latitude[combination$latitude,], timeExtent = as.POSIXct(allRanges$time[combination$time,], origin = '1970-01-01'), parentBrick = parentBrick, observations = parentBrick$observations)
  })

  if (gridCreationFlag) {
    parentBrick$M <- 1
    parentBrick$logLik <- NULL
    class(parentBrick) <- "Spacetimegrid"
    return(parentBrick)
  }
  invisible()
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
#' @param gridObj Spacetimegrid object
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

.addLayer <- function(gridObj, latBreaks, lonBreaks, timeBreaks) {
  gridObj$M <- gridObj$M+1
  newBreaks <- list(longitude = sort(c(tail(gridObj$breaks, n = 1)[[1]]$longitude, lonBreaks)), latitude = sort(c(tail(gridObj$breaks, n = 1)[[1]]$latitude, latBreaks)), time = sort(c(tail(gridObj$breaks, n = 1)[[1]]$time, timeBreaks)))
  gridObj$breaks <- c(gridObj$breaks, list(newBreaks))
  tipAddresses <- .tipAddresses(gridObj)
  lapply(tipAddresses, FUN = function(tipAddress) {
    .SpacetimegridConstructor(parentBrick = tipAddress, lonBreaks = lonBreaks, latBreaks = latBreaks, timeBreaks = timeBreaks)
    NULL
  })
  invisible()
}

#' Add breaks in space-time grid.
#'
#' This function is used to refine a resolution grid.
#'
#' @param gridObj Spacetimegrid object
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

# .sameGridSection <- function(x, y, brickObj) {
#   testResultTime <- identical(match(TRUE, brickObj$breaks[[m]][["time"]] > .getSpacetimeDim(x, "time")), match(TRUE, brickObj$breaks[[m]][["time"]] > .getSpacetimeDim(y, "time")))
#   if (!testResultTime) {
#     return(FALSE)
#   }
#   testResultsSpace <- sapply(c("longitude", "latitude"), FUN = function(dimensionName) {
#     xUpperPos <- match(TRUE, brickObj$breaks[[dimensionName]] > .getSpacetimeDim(x, dimensionName))
#     yUpperPos <- match(TRUE, brickObj$breaks[[dimensionName]] > .getSpacetimeDim(y, dimensionName))
#     identical(xUpperPos, yUpperPos)
#   })
#   all(testResultsSpace)
# }

.sameGridSection <- function(x, y, brickObj) {
  result <- TRUE
  if (!all(sapply(list(x,y), FUN = .coorWithinBrick, brickObj = brickObj))) {
    result <- FALSE
  }
  result
}

.coorWithinBrick <- function(spacetimeObj, brickObj) {
  timeTest <- (index(spacetimeObj)[[1]] >= min(brickObj$dimensions$time)) & (index(spacetimeObj)[[1]] < max(brickObj$dimensions$time))
  if (!timeTest) {
    return(FALSE)
  }
  dimNames <- c("longitude", "latitude")
  coordVec <- spacetimeObj@sp@coords[1, ]
  names(coordVec) <- dimNames
  spaceTests <- sapply(dimNames, function(dimName) (coordVec[[dimName]] >= min(brickObj$dimensions[[dimName]])) & (coordVec[[dimName]] < max(brickObj$dimensions[[dimName]])))
  if (!all(spaceTests)) {
    return(FALSE)
  }
  TRUE
}

.SpacetimebrickConstructor <- function(lonExtent, latExtent, timeExtent, parentBrick = NULL, observations = NULL) {
  if (is.null(parentBrick)) {
    parentBrick <- emptyenv()
    currentDepth <- 0
  }
  else {
    currentDepth <- parentBrick$depth + 1
  }

  if (!is.null(observations)) {
    observationsInRegion <- subset(observations, lonExtent = lonExtent, latExtent = latExtent, timeExtent = timeExtent)
  }
  brickEnvironment <- new.env(parent = parentBrick)
  brickEnvironment$depth <- currentDepth
  brickEnvironment$dimensions <- list(longitude = lonExtent, latitude = latExtent, time = timeExtent)
  brickEnvironment$observations <- observationsInRegion
  brickEnvironment$knotPositions <- NULL
  brickEnvironment$childBricks <- NULL

  if (!identical(parentBrick, emptyenv())) {
    childAddresses <- parentBrick$childBricks
    incAddresses <- c(childAddresses, brickEnvironment)
    parentBrick$childBricks <- incAddresses
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

.getM <- function(gridObj) {
  length(gridObj$breaks)
}

subset.STI <- function(x, latExtent, lonExtent, timeExtent) {
  timeIndices <- (index(x) <= max(timeExtent)) & (index(x) > min(timeExtent))
  lonIndices <- (x@sp@coords[, 1] <= max(lonExtent)) & (x@sp@coords[, 1] > min(lonExtent))
  latIndices <- (x@sp@coords[, 2] <= max(latExtent)) & (x@sp@coords[, 2] > min(latExtent))
  x[which(timeIndices & lonIndices & latIndices)]
}

.tipAddresses <- function(spacetimebrickObj) {
  if (!is.null(spacetimebrickObj$childBricks)) {
    return(unlist(lapply(spacetimebrickObj$childBricks, FUN = .tipAddresses)))
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
#' @param gridObj Spacetimegrid object
#' @param knotsList list of knot positions in spacetime format. If it has length 1, it will be assumed to be the same for all dimensions
#' @param r scalar or vector indicating the number of knots in each section of the grid at each resolution. Ignored if knotsList is specified. If it is a scalar, it is assumed that the number of knots is the same for each section of the grid at every resolution.
#' @param tuningPara if knots are placed automatically, they will be placed on the edges of a rectangle within each section. That rectangle is obtained by scaling the edges of each section by a coefficient equal to (1-2*tuningPara).
#'
#' @details This function operates with side-effects! It modifies gridObj, specifying the knotPositions component of each Spacetimebrick object. There might be advantages to having knots close to the edges, hence the automatic scaling scheme.
#'
#' @return Nothing: it directly modifies gridObj. In other words, it relies on side-effects.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

.addKnots <- function(gridObj, knotsList, argsForRandomKnots) {
  if (("STI" %in% class(knotsList))) { # If knotsList is NULL, class(knotsList) will return "NULL" (the class of NULL if "NULL", not inexistent.)
    knotsList <- replicate(.getM(gridObj)+1, expr = knotsList, simplify = FALSE) # In this situation, the knots are the same at every resolution. Might not be the best idea...
  }
  if (is.list(knotsList) & (length(knotsList) == 1)) {
    knotsList <- replicate(.getM(gridObj)+1, expr = knotsList[[1]], simplify = FALSE)
  }
  if (!is.null(knotsList)) {
    .populateKnots(gridObj, knotsList)
    invisible()
  }
  if (is.null(argsForRandomKnots)) {
    argsForRandomKnots <- list()
  }
  argsForRandomKnots$gridObj <- gridObj
  do.call(".populateKnotsRandomGrid", args = argsForRandomKnots)
  invisible()
}

.populateKnots <- function(gridObj, knotsList) {
  gridObj$knotPositions <- knotsList[[1]]
  if (!is.null(gridObj$childBricks)) {
    lapply(gridObj$childBricks, .populateKnots, knotsList = knotsList[-1])
  }
  invisible()
}

.setPointsOnRectangle <- function(numPoints, upperLeftCoor, width, height) {
  returnCoord <- function(distFromOrigin, upperLeftCoor, width, height) {
    upperRightCoor <- upperLeftCoor + c(width, 0)
    lowerRightCoor <- upperRightCoor - c(0, height)
    lowerLeftCoor <- upperLeftCoor - c(0, height)

    if (distFromOrigin < width) {
      return(upperLeftCoor + c(distFromOrigin,0))
    }
    if ((distFromOrigin >= width) & (distFromOrigin < width+height)) {
      return(upperRightCoor - c(0, distFromOrigin-width))
    }
    if ((distFromOrigin >= width+height) & (distFromOrigin < 2*width+height)) {
      return(lowerRightCoor - c(distFromOrigin-width-height, 0))
    }
    if ((distFromOrigin >= 2*width+height) & (distFromOrigin < 2*width+2*height)) {
      return(lowerLeftCoor + c(0, distFromOrigin - 2*width-height))
    }
    stop("Error: distance from origin is out of bounds! \n")
  }
  perimeter <- 2*(width+height)
  distBetweenPoints <- perimeter/(numPoints+1)
  distances <- seq(from = distBetweenPoints/2, to = perimeter-distBetweenPoints/2, length.out = numPoints)
  SpatialPoints(t(sapply(distances, returnCoord, upperLeftCoor = upperLeftCoor, width = width, height = height)))
}

.getAllParentAddresses <- function(brickObj) {

  brickObjList <- vector("list", length = brickObj$depth + 1)
  brickObjList[[brickObj$depth + 1]] <- brickObj

  if (!(brickObj$depth == 0)) {
    brickObjList <- .recurseParent(parent.env(brickObj), brickObjList)
  }
  brickObjList
}

.recurseParent <- function(currentBrick, brickObjList) {

  brickObjList[[currentBrick$depth + 1]] <- currentBrick
  if (!(currentBrick$depth == 0)) {
    brickObjList <- .recurseParent(parent.env(currentBrick), brickObjList)
  }
  brickObjList
}

logLik.Spacetimegrid <- function(gridObj) {
  gridObj$logLik
}

# To make computations lighter, we assume that knot locations at the finest resolution match observations.
# It's also important that no two knot locations match.
# The number of knots should be highest at the finest resolution.
# Knots are placed solely within regions where observations are found.
# The basic scheme will place the knots randomly in the spatial coordinates, drawing the coordinates from a uniform distribution. Since time is discrete, the time coordinate of each knot will also be drawn at random, but from a discrete distribution.
# The number of knots will be determined with a scaling constant, e.g. for M=3, (1/4, 1/2, 3/4, 1)*n. Any strictly positive function strictly increasing defined over 0 to M would be valid, with the additional condition that it takes value 1 at M, e.g. (1/M+1)*(x+1), sqrt((x+1)/(M+1)), log((x+1)/(M+1))+1, exp(theta(x - M))

.populateKnotsRandomGrid <- function(gridObj, numKnotsFun = NULL, seed = NULL) {
  set.seed(seed)
  if (is.null(numKnotsFun)) {
    numKnotsFun <- function(m) 1/(gridObj$M+1)*(m+1)
  }
  if (is.null(seed)) {
    seed <- 42
  }
  if (gridObj$M == 0) {
    gridObj$knotPositions <- gridObj$observations$sp@coords
    invisible()
  }
  numKnots <- ceiling(numKnotsFun(0)*length(gridObj$observations))
  if (numKnots > 0) {
    randomizeKnotsInBrick(gridObj, numKnots)
    lapply(gridObj$childBricks, .populateKnotsRandomBrick, numKnotsFun = numKnotsFun, M = gridObj$M)
  }
  invisible()
}

.populateKnotsRandomBrick <- function(brickObj, numKnotsFun, M) {
  if (brickObj$depth == M) {
    brickObj$knotPositions <- brickObj$observations@sp@coords
  }
  numKnots <- ceiling(numKnotsFun(brickObj$depth)*length(brickObj$observations))
  if (numKnots > 0) {
    randomizeKnotsInBrick(brickObj, numKnots)
    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, .populateKnotsRandomBrick, numKnotsFun = numKnotsFun, M = M)
    }
  }
  invisible()
}

randomizeKnotsInBrick <- function(brickObj, numKnots) {
  shortRunif <- function(n, bounds) runif(n, bounds[1], bounds[2])
  brickKnotsLon <- shortRunif(numKnots, brickObj$dimensions$longitude)
  brickKnotsLat <- shortRunif(numKnots, brickObj$dimensions$latitude)
  brickKnotsTimeIndices <- sample(1:length(brickObj$observations@time), numKnots, replace = TRUE)
  knotPositions <- STI(sp = SpatialPoints(cbind(brickKnotsLon, brickKnotsLat)), time = brickObj$observations@time[brickKnotsTimeIndices])
  brickObj$knotPositions <- knotPositions
  invisible()
}
