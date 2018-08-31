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

MRA_INLA <- function(data, covFct, gridObj, numKnots, hyperpriorFunList) {
 # TO_DO
}

.logLikFun <- function() {
  # TO_DO
}

.logGMRFprior <- function(paraValuesList, SigmaValuesList) {
  sum(mapply(paraValues = paraValuesList, varPar = SigmaValuesList, FUN = function(paraValues, varPar) {
    if (!is.matrix(varPar)) {
      varPar <- varPar*diag(length(paraValues))
    }
    mvtnorm::dmvnorm(x = paraValues, sigma = varPar, log = TRUE)
  }))
}

setupGrid <- function(lonNewBreaksList, latNewBreaksList, timeNewBreaksList, observations = NULL, knotsList = NULL, r = rep(10, length(lonNewBreaksList)+1), covFct, ...) {
  gridForMRA <- .SpacetimegridConstructor(parentBrick = NULL, lonBreaks = lonNewBreaksList[[1]], latBreaks = latNewBreaksList[[1]], timeBreaks = timeNewBreaksList[[1]], observations = observations)
  if (length(lonNewBreaksList) == 1) {
    return(gridForMRA)
  }
  lapply(seq_along(lonNewBreaksList)[-1], FUN = function(resolutionIndex) {
    .addLayer(gridForMRA, latBreaks = latNewBreaksList[[resolutionIndex]], lonBreaks = lonNewBreaksList[[resolutionIndex]], timeBreaks = timeNewBreaksList[[resolutionIndex]])
  })
  .addKnots(gridForMRA, r = r, ...)
  .getKmatsAndLocalFunctions(gridForMRA, covFct)
  gridForMRA
}

# The next function returns the K matrices, b and v functions for all zones of the pre-defined grids.

.getKmatsAndLocalFunctions <- function(gridObj, covFct) {
  gridObj$K <-  Matrix::chol2inv(Matrix::chol(covFct(gridObj$knotPositions, gridObj$knotPositions)))
  gridObj$vFun <- covFct
  gridObj$bFun <- .getBfun(gridObj)

  if (!is.null(gridObj$childBricks)) {
    lapply(gridObj$childBricks, FUN = .genBrickFcts)
  }
  invisible()
}

.genBrickFcts <- function(brickObj) {
  .updateVfun(brickObj)
  .getBfun(brickObj)
  .getKmat(brickObj)
  if (!is.null(brickObj$childBricks)) {
    lapply(brickObj$childBricks, .genBrickFcts)
  }
  invisible()
}

.updateVfun <- function(brickObj) {
  brickObj$vFun <- function(spacetime1, spacetime2) {
    if(!.sameGridSection(spacetime1, spacetime2, parent.env(brickObj))) {
      return(0)
    }
    parent.env(brickObj)$vFun(spacetime1, spacetime2) - t(parent.env(brickObj)$bFun(spacetime1))%*%parent.env(brickObj)$K%*%parent.env(brickObj)$bFun(spacetime2)
  }
  invisible()
}

.getBfun <- function(brickObj) {
  brickObj$bFun <- function(spacetimeCoord) {
    if (!.coorWithinBrick(spacetimeCoord, brickObj)) {
      return(rep(0, Npoints(brickObj$knotPositions)))
    }
    sapply(1:Npoints(brickObj$knotPositions), FUN = function(knotIndex) {
      brickObj$vFun(spacetimeCoord, brickObj$knotPositions[knotIndex])
    })
  }
}

.getKmat <- function(brickObj) {
  numKnots <- Npoints(brickObj$knotPositions)
  KmatrixFinal <- matrix(0, numKnots, numKnots)
  # Handling the off diagonal elements
  diagMat <- diag(numKnots)
  mapply(rowIndex = row(diagMat)[lower.tri(diagMat)], colIndex = col(diagMat)[lower.tri(diagMat)], FUN = function(rowIndex, colIndex) {
    KmatrixFinal[rowIndex, colIndex] <<- brickObj$vFun(brickObj$knotPositions[rowIndex], brickObj$knotPositions[colIndex])
    invisible()
  }, SIMPLIFY = FALSE)
  KmatrixFinal <- KmatrixFinal + t(KmatrixFinal)
  diag(KmatrixFinal) <- sapply(1:numKnots, FUN = function(x) brickObj$vFun(brickObj$knotPositions[x], brickObj$knotPositions[x]))
  brickObj$K <- KmatrixFinal
  invisible()
}

spacetimeListConvertToPoints <- function(valuesList, timeValues = NULL, regular = FALSE) {
  if (regular) {
    if (class(valuesList[[1]]) == "SpatialPoints") {
      return(spacetime::STF(sp = valuesList[[1]], time = timeValues))
    }
    valuesFrames <- lapply(valuesList, FUN = function(x) x@data)
    dataReformat <- do.call("rbind", valuesFrames)
    return(spacetime::STFDF(valuesList[[1]], time = timeValues, data = dataReformat))
  }
  allPoints <- do.call(what = raster::bind, valuesList)
  lengthVec <- sapply(valuesList, length)
  extendedTimeValues <- rep(timeValues, lengthVec)

  if (class(allPoints) == "SpatialPoints") {
    return(spacetime::STI(sp = allPoints, time = extendedTimeValues))
  }
  spacetime::STIDF(sp = allPoints, time = extendedTimeValues, data = allPoint@data)
}

Npoints <- function(spacetimeObj) {
  nrow(spacetimeObj@sp@coords)
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

