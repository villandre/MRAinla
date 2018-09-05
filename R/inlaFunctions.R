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
  .computeWmats(gridObj = gridForMRA, covFct = covFct)
  .setBtips(gridForMRA)
  .setAtildeTips(gridForMRA)
  .recurseA(gridForMRA)
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

.computeWmats <- function(covFct, gridObj) {
  gridObj$Kinverse <- covFct(gridObj$knotPositions, gridObj$knotPositions)
  gridObj$K <- Matrix::chol2inv(Matrix::chol(gridObj$Kinverse))
  gridObj$W <- list(gridObj$K)
  if (!is.null(gridObj$childBrick)) {
    lapply(gridObj$childBricks, .computeWchildBrick, covFct = covFct)
  }
  invisible()
}

.computeWchildBrick <- function(brickObj, covFct) {
  brickObj$W <- vector("list", brickObj$depth+1) # Depth is from 0 to M.
  brickPointers <- .getAllParentAddresses(brickObj) # brickPointers includes the address of brickObj from depth 0 at position 1 to depth brickObj$depth at position brickObj$depth+1.
  brickObj$W[[1]] <- covFct(brickObj$knotPositions, brickPointers[[1]]$knotPositions)

  for (l in (seq_along(brickPointers))) {
    brickPointer <- brickPointers[[l+1]]
    covEvalation <- covFct(brickObj$knotPositions, brickPointer$knotPositions)
    headWlist <- brickObj$W[1:l]
    Kmatrices <- lapply(brickPointers[1:l], function(aPointer) aPointer$K)
    tailWlist <- lapply(brickPointer$W[1:l], t)
    secondTerm <- Reduce("+", mapply(x = headWlist, y = Kmatrices, z = tailWlist, FUN = function(x,y,z) x%*%y%*%z, SIMPLIFY = FALSE))
    brickObj$W[[l+1]] <- covEvalation - secondTerm
  }
  brickObj$Kinverse <- tail(brickObj$W, n = 1)[[1]]
  brickObj$K <-  Matrix::chol2inv(Matrix::chol(brickObj$Kinverse))
  if (!is.null(brickObj$childBricks)) {
    lapply(brickObj$childBricks, .computeWchildBrick, covFct = covFct)
  }
  invisible()
}

.setDUtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  transferW <- function(tipBrick) {
    SigmaMat <- tail(tipBrick$W, n = 1)[[1]]
    tipBrick$SigmaMat <- SigmaMat
    tipBrick$d <- log(det(SigmaMat))
    tipBrick$u <- t(tipBrick$observations)%*%Matrix::chol2inv(Matrix::chol(SigmaMat))%*%tipBrick$observations
    tipBrick$KtildeInverse <- tipBrick$Kinverse + t(tipBrick$bMat)%*%tipBrick$Vinverse%*%tipBrick$bMat
    invisible()
  }
  lapply(allTips, transferW)
  invisible()
}

.setDUinside <- function(brickObj) {
  if (!is.null(gridObj$childBricks)) {
    if (is.null(gridObj$childBricks[[1]]$SigmaMat)) {
      lapply(gridObj$childBricks, .setDUinside)
    }
    dSecondPart <- -log(det(brickObj$Kinverse)) + sapply(brickObj$childBricks, function(x) x$d)
    dFirstPart <- log(det(brickObj$KtildeInverse))
    brickObj$d <- dFirstPart + dSecondPart
    uSecondPart <- sapply(brickObj$childBricks, function(x) x$u)
    # uFirstPart <- t(brickObj$omega)%*%brickObj$Ktilde%*%brickObj$omega
    uFirstPart <- .colSums(brickObj$omega * (brickObj$Ktilde %*% brickObj$omega), m = length(brickObj$omega), n = 1) # A faster way to compute quadratic forms.
    gridObj$SigmaMat <- gridObj$bMat%*%gridObj$K%*%t(gridObj$bMat) + gridObj$V
  }
  invisible()
}

.getKtildeInverse <- function(brickObj) {
  brickObj$Kinverse
}

.setBtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  lapply(allTips, FUN = function(x) {
    x$BmatList <- x$Wmat[1:(x$depth+1)]
    invisible()
  }) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.setAtildeTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)

  lapply(allTips, FUN = function(x) {
    x$Atildemm <- t(x$BmatList[[gridObj$M]]) %*% x$K %*% x$BmatList[[gridObj$M]] ## BmatList has length M+1, but element M+1 is NULL. gridDepth is equal to M, since there's a resolution 0. currentDepth in this situation is the second to last element, which was defined when .setBtips was called.
    invisible()
  }) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}
.recurseA <- function(brickObj) {
  if (is.null(brickObj$childBricks)) {
    return(NULL)
  }

  if (is.null(brickObj$childBricks[[1]]$Atildemm)) {
    lapply(brickObj$childBricks, .recurseA)
  }
  AtildeMatrices <- lapply(brickObj$childBricks, function(x) x$Atildemm)
  brickObj$Amm <- Reduce("+", AtildeMatrices)
  brickObj$KtildeInverse <- brickObj$Kinverse + brickObj$Amm
  brickObj$Atildemm <- brickObj$Amm - brickObj$Amm %*% brickObj$Ktilde %*% brickObj$Amm
  invisible()
}
