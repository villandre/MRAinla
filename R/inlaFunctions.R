#' Fit the MRA model for spatiotemporal data
#'
#' Implementation of the spatiotemporal of Villandre et al. 2018.
#'
#' @param spaceTimeList list of SpatialPointsDataFrame objects
#' @param spaceTimeCovFct function with at least two arguments, spacetimeCoord1 and spacetimeCoord2
#' @param M number of embedded resolutions
#' used in proposing transitions in the space of between-cluster phylogenies
#' @param gridRasterList list of rasters splitting the domain into embedded resolutions
#' @param numKnots either the number of knots in all elements of the grids, or a vector of size equal to
#' the length of gridRasterList giving the number of knots for each resolution, or a list of vectors giving
#' the number of knots in each region of the grid for each resolution
#' @param hyperPriorFunList named list of functions with one argument specifying the hyperprior distributions
#'
#' @details Nothing to say for now.
#'
#' @return A list with two components:
#' \itemize{
#'  \item{chain:} {list with each element itself a list giving the sampled parameter values.
#' betweenTransMatListIndex and withinTransMatListIndex correspond to the index of the assumed mean
#' branch length in meanBetweenBranchVec and meanWithinBranchVec, respectively}
#'  \item{MAPestimate:}{ vector giving the maximum posterior probability cluster membership indices estimate}
#' }
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

MRA_INLA <- function(spaceTimeList, spaceTimeCovFct, M, gridRasterList, numKnots, hyperpriorFunList) {

}

.logLikFun <- function(spaceTimeList, meanValueList, sigmaErrorValue) {
  sum(mapply(grid = spaceTimeList, meanValues = meanValueList, FUN = function(grid, meanValues) {
    sum(dnorm(x = grid@data, mean = meanValues, sd = sigmaErrorValue, log = TRUE))
  }))
}

.logGMRFprior <- function(paraValuesList, SigmaValuesList) {
  sum(mapply(paraValues = paraValuesList, varPar = SigmaValuesList, FUN = function(paraValues, varPar) {
    if (!is.matrix(varPar)) {
      varPar <- varPar*diag(length(paraValues))
    }
    mvtnorm::dmvnorm(x = paraValues, sigma = varPar, log = TRUE)
  }))
}

.setupVrecursionStep <- function(spacetimegridObj, v, baseVec1, baseVec2, K) {
  function(spaceTime1, spaceTime2) {
    if(!.sameGridSection(spaceTime1, spaceTime2, spacetimeGridObj)) {
      return(0)
    }
    v - baseVec1%*%K%*%baseVec2
  }
}

.initVrecursion <- function(knots, covFct) {
  function(spaceTime1, spaceTime2) {
    initBList <- list(covFct(spaceTime1, knots), covFct(spaceTime2, knots))
    initInvK <- covFct(knots, knots)
    initK <- Matrix::chol2inv(Matrix::chol(initInvK))
    initV <- covFct(spaceTime1, spaceTime2)
    list(v = initV, K = initK, bList = initBList)
  }
}

.setupFullrecursion <- function(gridList, knotsList, covFct) {
  initialVfun <- initVrecursion(knots = knotsList[[1]], covFct = covFct)
  function(spaceTime1, spaceTime2) {
    currentValues <- initialVfun(spaceTime1, spaceTime2)
    lapply(seq_along(gridList), FUN = function(resIndex) {
      incrementedVfun <- setupVrecursionStep(spacetimegridObj = gridList[[resIndex]], v = currentValues$v, baseVec1 = currentValues$bList[[1]], baseVec2 = currentValues$bList[[2]], K = currentValues$K)
      newV <- incrementedVfun(spaceTime1, spaceTime2)
      newBs <- lapply(c(spaceTime1, spaceTime2), FUN = getBvec, knotPositions = knotsInRegion, vFun = incrementedVfun)
      newK <- getKmat(knotPositions = knotsInRegion , vFun = incrementedVfun)
    })
  }
}

.getBvec <- function(spaceTime, knotPositions, vFun) {
  sapply(knotPositions, FUN = function(knotPos) {
    vFun(spaceTimeCoord, knotPos)
  })
}

.getKmat <- function(knotPositions, vFun) {
  KmatrixFinal <- matrix(0, nrow(knotPositions), nrow(knotPositions))
  # Handling the off diagonal elements
  mapply(rowIndex = row(diag(nrow(knotPositions)))[lower.tri(diag(nrow(knotPositions)))], colIndex = col(diag(nrow(knotPositions)))[lower.tri(diag(nrow(knotPositions)))], FUN =function(rowIndex, colIndex) {
    vFun(knotPositions[rowIndex], knotPositions[colIndex])
    invisible(NULL)
  }, SIMPLIFY = FALSE)
  diag(KmatrixFinal) <- sapply(1:nrow(knotPositions), FUN = function(x) {
    vFun(knotPositions[x], knotPositions[x])
  })
}

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
  class(grids) <- "Spacetimegrid"
  grids
}

splitGridSection <- function(gridObject, dimension = c("lon", "lat", "time"), breakPos) {
  gridObject[[dimension]]$breaks <- sort(c(gridObject[[dimension]]$breaks, breakPos))
  gridObject
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
  cat("Longitude extent: ", x$longitude$extent , "\n")
  cat("Latitude extent: ", x$latitude$extent, "\n")
  cat("Time extent: ", x$time$extent, "\n \n")
  cat("Grid dimensions: (", length(x$longitude$breaks) + 1, " ", length(x$latitude$breaks)+1, " ", length(x$time$breaks) + 1, ")\n")
  cat("Note: This gives the number of regions, which may have different volumes.\n\n")
  invisible(NULL)
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
    combinedSP <- bind(combinedSP, knotSP)
    plotColours <- combinedSP@data$identity
  }
  plot(combinedSP, col = as.character(plotColours), xlab = "Longitude", ylab = "Latitude", xlim = x$longitude$extent, ylim = x$latitude$extent, cex = 0.8, pch = 18)
  abline(v = c(x$longitude$extent[1], x$longitude$breaks, x$longitude$extent[2]), lty = 3, col = "blue4")
  abline(h = c(x$latitude$extent[1], x$latitude$breaks, x$latitude$extent[2]), lty = 3, col = "blue4")
}

# observations are coded as SpatialPointsDataFrame with 3D coordinates, longitude-latitude-time

getPointsInRegion <- function(sptGrid, pointSpacetimeCoor, observations) {
  mapply(c("longitude", "latitude", "time"), coorPos = 1:3, FUN = function(dimensionName, coorPos) {
    upperPos <- match(TRUE, sptGrid[[dimensionName]]$breaks > pointSpacetimeCoor@coords[1, coorPos])
    coorRange <- c(sptGrid[[dimensionName]]$breaks[upperPos - 1], sptGrid[[dimensionName]]$breaks[upperPos])
    observations <<- observations[(observations@coords[ , coorPos] > coorRange[[1]]) & (observations@coords[ , coorPos] < coorRange[[2]])]
    invisible(NULL)
  }, SIMPLIFY = FALSE)
  observations
}

.sameGridSection <- function(x, y, spacetimegridObj) {
  sameSection <- TRUE
  testResults <- mapply(c("longitude", "latitude", "time"), coorPos = 1:3, FUN = function(dimensionName, coorPos) {
    xUpperPos <- match(TRUE, spacetimeGridObj[[dimensionName]]$breaks > x@coords[1, coorPos])
    yUpperPos <- match(TRUE, spacetimeGridObj[[dimensionName]]$breaks > y@coords[1, coorPos])
    identical(xUpperPos, yUpperPos)
  })
  if (!all(testResults)) {
    sameSection <- FALSE
  }
  sameSection
}
