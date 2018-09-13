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

.logGMRFprior <- function(paraValuesList, SigmaValuesList) {
  sum(mapply(paraValues = paraValuesList, varPar = SigmaValuesList, FUN = function(paraValues, varPar) {
    if (!is.matrix(varPar)) {
      varPar <- varPar*diag(length(paraValues))
    }
    mvtnorm::dmvnorm(x = paraValues, sigma = varPar, log = TRUE)
  }))
}

#' Set up a multi-resolution spatiotemporal grid for MRA estimation
#'
#' The spatiotemporal grid has one layer per resolution (including one for resolution 0).
#'
#' @param lonNewBreaksList list of longitude breaks, one element per resolution, breaks are specified incrementally (see Details)
#' @param latNewBreaksList list of latitude breaks, one element per resolution, breaks are specified incrementally (see Details)
#' @param timeNewBreaksList list of time breaks, one element per resolution, breaks are specified incrementally (see Details)
#' @param observations a Spacetime object giving the position and values of the observed data
#' @param knotsList (optional) list of Spacetime objects giving knot positions at each resolution. If left empty, knots are placed automatically
#' @param r (optional) If knots are placed automatically, a vector indicating the number of knots to create in each region
#' @param covFct a function with two arguments, each of which a Spacetime object with one point
#' @param ... Additional parameters for the automated knot placement scheme, see .populateKnots
#'
#' @details Since resolutions are nested, breaks should be specified incrementally, i.e. list(c(1, 2, 3, 4), c(1.5, 2.5, 3.5)) will create three grids, one for resolution 0 with no breaks, but coordinates ranging from 1 to 4, one for resolution 1 with additional breaks 2 and 3, and one for resolution 2 with additional breaks 1.5, 2.5, and 3.5.
#' The knot placement scheme places knots close to the boundaries of each spacetime brick, with tuningPara (if specified) determining how close in terms of the proportion of the range knots should be placed to spatial boundaries. For time, knots are placed at the extremities, keeping in mind that regions are closed on the left and open on the right, e.g \[2008-01-01, 2009-02-03).
#'
#' @return A Spacetimegrid object.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export


setupGrid <- function(lonNewBreaksList, latNewBreaksList, timeNewBreaksList, observations = NULL, knotsList = NULL, r = rep(10, length(lonNewBreaksList)+1), covFct, ...) {
  gridForMRA <- .SpacetimegridConstructor(parentBrick = NULL, lonBreaks = lonNewBreaksList[[1]], latBreaks = latNewBreaksList[[1]], timeBreaks = timeNewBreaksList[[1]], observations = observations)
  if (length(lonNewBreaksList) == 1) {
    return(gridForMRA)
  }
  lapply(seq_along(lonNewBreaksList)[-1], FUN = function(resolutionIndex) {
    .addLayer(gridForMRA, latBreaks = latNewBreaksList[[resolutionIndex]], lonBreaks = lonNewBreaksList[[resolutionIndex]], timeBreaks = timeNewBreaksList[[resolutionIndex]])
  })
  .addKnots(gridForMRA, r = r, ...)

  .computeWmats(gridObj = gridForMRA, covFct = covFct)
  .setBtips(gridForMRA)
  .setSigmaTips(gridForMRA)
  .setAtildeTips(gridForMRA)
  .recurseA(gridForMRA)
  .setOmegaTildeTips(gridForMRA)
  .recurseOmega(gridForMRA)

  .setDtips(gridForMRA)
  .recurseD(gridForMRA)

  .setUtips(gridForMRA)
  .recurseU(gridForMRA)

  gridForMRA$logLik <- gridForMRA$d + gridForMRA$u
  gridForMRA
}

spacetimeListConvertToPoints <- function(valuesList, timeValues = NULL) {
  allPoints <- do.call(what = raster::bind, valuesList)
  lengthVec <- sapply(valuesList, length)
  extendedTimeValues <- rep(timeValues, lengthVec)

  if (class(allPoints) == "SpatialPoints") {
    return(spacetime::STI(sp = allPoints, time = extendedTimeValues))
  }
  spacetime::STIDF(sp = SpatialPoints(allPoints@coords), time = extendedTimeValues, data = allPoints@data) # Having the sp argument being a SpatialPointsDataFrame object produces an error. Stripping it of its data component then re-specifying it fixes the issue.
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

.computeWmats <- function(gridObj, covFct) {
  gridObj$Kinverse <- covFct(gridObj$knotPositions, gridObj$knotPositions)
  gridObj$K <- Matrix::chol2inv(Matrix::chol(gridObj$Kinverse))
  gridObj$WmatList <- list(gridObj$K)
  if (!is.null(gridObj$childBricks)) {
    lapply(gridObj$childBricks, .computeWchildBrick, covFct = covFct)
  }
  invisible()
}

.computeWchildBrick <- function(brickObj, covFct) {
  brickObj$WmatList <- vector("list", brickObj$depth+1) # Depth is from 0 to M.
  brickPointers <- .getAllParentAddresses(brickObj) # brickPointers includes the address of bricks from depth 0 (position 1) to depth brickObj$depth (position brickObj$depth+1).
  brickObj$WmatList[[1]] <- covFct(brickObj$knotPositions, brickPointers[[1]]$knotPositions)
  m <- brickObj$depth

  for (l in 0:m) { # First element in brickPointers is the top brick, i.e. at resolution 0.
    brickPointer <- brickPointers[[l+1]]
    covEvaluation <- covFct(brickObj$knotPositions, brickPointer$knotPositions)
    secondTerm <- 0
    if (l > 0) {
      headWlist <- brickObj$WmatList[(0:l-1) + 1]
      Kmatrices <- lapply(brickPointers[(0:l-1) + 1], function(aPointer) aPointer$K)
      tailWlist <- lapply(brickPointer$WmatList[(0:l-1) + 1], t)

      secondTerm <- Reduce("+", mapply(x = headWlist, y = Kmatrices, z = tailWlist, FUN = function(x,y,z) x %*% y %*% z, SIMPLIFY = FALSE))
    }

    brickObj$WmatList[[l+1]] <- covEvaluation - secondTerm
  }
  if (brickObj$depth < .getTopEnvirAddress(brickObj)$M) {
    brickObj$Kinverse <- tail(brickObj$WmatList, n = 1)[[1]]
    brickObj$K <- Matrix::chol2inv(Matrix::chol(brickObj$Kinverse))
  }
  if (!is.null(brickObj$childBricks)) {
    lapply(brickObj$childBricks, .computeWchildBrick, covFct = covFct)
  }
  invisible()
}

.setBtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  lapply(allTips, FUN = function(x) {
    x$BmatList <- x$WmatList[1:(x$depth+1)]
    invisible()
  }) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.setSigmaTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  lapply(allTips, FUN = function(x) {
    x$Sigma <- x$WmatList[[gridObj$M+1]]
    x$SigmaInverse <- solve(x$Sigma)
    invisible()
  }) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.setAtildeTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  m <- gridObj$M

  modifyTip <- function(x) {

    x$Atilde <- vector('list', x$depth+1) # First level of the hierarchy is for k in Atilde^{k,l}.
    x$Atilde <- lapply(0:x$depth, function(k) vector('list', x$depth)) # Second level is for l in Atilde^{k,l}, same order as for k.
    indexGrid <- expand.grid(k = 0:(m-1), l = 0:(m-1)) # We don't need A^{M,M}_{j_1, ..., j_M}
    indexGrid <- subset(indexGrid, subset = k >= l)

    for (i in 1:nrow(indexGrid)) {
      k <- indexGrid[i, 1]
      l <- indexGrid[i, 2]
      x$Atilde[[k + 1]][[l + 1]] <- t(x$BmatList[[k + 1]]) %*% x$SigmaInverse %*% x$BmatList[[l + 1]] ## BmatList has length M+1, but element M+1 is NULL. gridDepth is equal to M, since there's a resolution 0. Element gridObj$M in this situation is the second to last element, which was defined when .setBtips was called.
      x$Atilde[[l + 1]][[k + 1]] <- t(x$Atilde[[k + 1]][[l + 1]])
    }
    invisible()
  }
  lapply(allTips, FUN = modifyTip) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.recurseA <- function(brickObj) {
  if (is.null(brickObj$childBricks)) {
    return(NULL)
  }

  if (is.null(brickObj$childBricks[[1]]$Atilde)) {
    lapply(brickObj$childBricks, .recurseA)
  }

  m <- brickObj$depth

  indexGrid <- expand.grid(k = 0:m, l = 0:m)
  indexGrid <- subset(indexGrid, subset = k >= l)

  Amatrices <- mapply(k = indexGrid$k, l = indexGrid$l, FUN = function(k,l) { # +1 to adjust indices for R.
    Reduce("+", lapply(brickObj$childBricks, function(x) {
      x$Atilde[[k + 1]][[l + 1]]
    }))
  }, SIMPLIFY = FALSE)
  brickObj$A <- .reformatList(matrixList = Amatrices, indexGrid = indexGrid)

  brickObj$KtildeInverse <- brickObj$Kinverse + brickObj$A[[m + 1]][[m + 1]] # +1 to adjust indices for R.
  brickObj$Ktilde <- Matrix::chol2inv(chol(brickObj$KtildeInverse))

  AtildeMatrices <- mapply(k = indexGrid$k, l = indexGrid$l, FUN = function(k, l) {
    brickObj$A[[k + 1]][[l + 1]] - brickObj$A[[k + 1]][[m + 1]] %*% brickObj$Ktilde %*% brickObj$A[[m + 1]][[l + 1]] # m+1 to adjust indices for R.
  }, SIMPLIFY = FALSE)
  brickObj$Atilde <- .reformatList(matrixList = AtildeMatrices, indexGrid = indexGrid)

  invisible()
}

.computeDtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)

  modifyTip <- function(x) {
    x$d <- log(det(x$Sigma))
    invisible()
  }
  lapply(allTips, FUN = modifyTip)
  invisible()
}

.recurseD <- function(brickObj) {
  if (is.null(brickObj$childBricks)) {
    return(invisible())
  }

  if (is.null(brickObj$childBricks[[1]]$d)) {
    lapply(brickObj$childBricks, .computeD)
  }

  brickObj$d <- log(det(brickObj$KtildeInverse)) - log(det(brickObj$Kinverse)) + sum(vapply(brickObj$childBricks, FUN.VALUE = 'numeric', function(childBrick) childBrick$d))
  invisible()
}

.setOmegaTildeTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  m <- gridObj$M

  modifyTip <- function(x) {

    x$omegaTilde <- replicate(n = x$depth+1, expr = vector('list', x$depth + 1), simplify = FALSE)

    for (k in 0:(gridObj$M-1)) {

      x$omegaTilde[[k+1]] <- t(x$BmatList[[k+1]]) %*% x$SigmaInverse %*% x$observations@data[ , 1]
    }
    invisible()
  }
  lapply(allTips, FUN = modifyTip) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.recurseOmega <- function(brickObj) {
  if (is.null(brickObj$childBricks)) {
    return(NULL)
  }

  if (is.null(brickObj$childBricks[[1]]$omegaTilde)) {
    lapply(brickObj$childBricks, .recurseOmega)
  }

  m <- brickObj$depth

  omegaMatrices <- lapply(0:m, FUN = function(k) {
    Reduce("+", lapply(brickObj$childBricks, function(x) {
      x$omegaTilde[[k+1]]
    }))
  }, SIMPLIFY = FALSE)

  brickObj$omegaTilde <- lapply(0:m, FUN = function(k) {
    brickObj$omega[[k+1]] - brickObj$A[[k+1]][[m+1]] %*% brickObj$Ktilde %*% brickObj$omega[[m+1]]
  })
  invisible()
}

.computeUtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  m <- gridObj$M

  modifyTip <- function(x) {
    x$u <- t(x$observations@data[ , 1]) %*% x$SigmaInverse %*% x$observations@data[ , 1] ## QUADRATIC FORM: MAY BE OPTIMIZED.
    invisible()
  }
  lapply(allTips, FUN = modifyTip)
  invisible()
}

.recurseU <- function(gridObj) {
  if (is.null(brickObj$childBricks)) {
    return(invisible())
  }

  if (is.null(brickObj$childBricks[[1]]$d)) {
    lapply(brickObj$childBricks, .computeD)
  }

  brickObj$u <- -t(brickObj$omega[[brickObj$depth+1]]) %*% brickObj$Ktilde %*% brickObj$omega[[brickObj$depth + 1]] + sum(vapply(brickObj$childBricks, FUN.VALUE = 'numeric', function(childBrick) childBrick$u)) ## QUADRATIC FORM, MIGHT BE POSSIBLE TO OPTIMIZE.
  invisible()
}

.reformatList <- function(matrixList, indexGrid) {
  m <- max(indexGrid[ , 1])
  newList <- replicate(m + 1, expr = vector('list', m + 1), simplify = FALSE)
  for (i in seq_along(matrixList)) {
    k <- indexGrid[i, 1]
    l <- indexGrid[i, 2]
    newList[[k + 1]][[l + 1]] <- matrixList[[i]]
    newList[[l + 1]][[k + 1]] <- t(newList[[k + 1]][[l + 1]])
  }
  newList
}
