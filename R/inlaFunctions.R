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

MRA_INLA <- function(covarianceFct, gridObj, numKnots, hyperpriorFunList) {
  hyperpriorNames <- names(hyperpriorFunList)
  covarianceParas <- hyperpriorNames[!(hyperpriorNames == "MEvar")]

  respAndCovValues <- gridObj$dataset@data
  psiFun <- function(psi) { # Values of psi must be named! The element for the measurement error variance must be named "MEvar".
    prod(sapply(seq_along(hyperpriorFunList), function(hyperIndex) hyperpriorFunList[hyperIndex](psi[hyperIndex])))
  }
  dataLikFun <- function(x, mean, variance) {
    dnorm(x = x, mean = mean, sd = sqrt(variance))
  }
  MRAlikFun <- function(gridObj, spacetimeTheta, hyperpara) {
    covFct <- function(spacetime1, spacetime2) {
      covarianceHyperValuesAsList <- as.list(hyperpara[covarianceParas])
      do.call("covarianceFct", args = c(list(spacetime1 = spacetime1, spacetime2 = spacetime2), covarianceHyperValuesAsList))
    }
    .populateObservations(gridObj, spacetimeTheta) # Maybe computational gains could be obtained if the thetas could be directly assigned to the different regions.
    computeLogLik(gridObj = gridObj, covFct = covarianceFct)
  }
  .makeMarginalPsiFun(gridObj, psiPriorFun = psiFun, MRAlikFun = MRAlikFun, dataLikFun = dataLikFun)
}

# MRAprecision will have to be coded as block diagonal to ensure tractability.

.makeMarginalPsiFun <- function(gridObj, psiPriorFun, MRAlikFun, dataLikFun) {
  function(psi) {
    MRAlikFun(gridObj = gridObj, spacetimeTheta = modeAndCurvature$mode, hyperpara = psi)
    modeAndCurvature <- .findModeAndCurvature(y = spacetimeData@data[, 1], covValues = spacetimeData@data[, -1], MEvar = psi$MEvar, MRAprecision = gridObj$SigmaInverse) # The MRA precision matrix is assumed to be SigmaInverse for the spatiotemporal brick at M=0.
    psiPriorFun(psi)*exp(gridObj$logLik)*dataLikFun(y, modeAndCurvature$mode, psi$MEvar)/mvnfast::dmvn(X = modeAndCurvature$mode, mu = modeAndCurvature$mode, sigma = modeAndCurvature$curvature, ncores = 1) # THE NUMBER OF CORES COULD BE CHANGED EVENTUALLY BUT LET'S SEE IF PERFORMANCE IN THIS STEP IS A REAL ISSUE BEFORE GOING MULTITHREADED. That'll probably result in a computational zero.
  }
}

.findModeAndCurvature <- function(y, covValues, startingValuesForTheta = rep(0, ncol(covValues) + 1 + nrow(covValues)), MEvar, MRAprecision) {
  InvFisherMatrix <- Matrix::chol2inv(chol(.FisherInfo(covValues = covValues, MEvar = MEvar, MRAprecision = MRAprecision)))
  shortScoreFun <- function(theta) {
    fixed <- theta[1:(nrow(covValues)+1)]
    ran <- theta[(nrow(covValues)+2):length(y)]
    .scoreFun(y = y, fixedCoefs = fixed, MEvar = MEvar, covValues = covValues, ranEffects = ran, MRAprecision = MRAprecision)
  }
  modeValue <- .FisherScoringAlgo(invFisherInfo = InvFisherMatrix, scoreFun = shortScoreFun)
  curvature <- .curvatureFun(covValues = covValues, MEvar = MEvar, MRAprecision = MRAprecision) # The curvature does not depend on the mode itself.
  list(mode = modeValue, curvature = curvature)
}

.scoreFun <- function(y, fixedCoefs, MEvar, covValues, ranEffects, MRAprecision) {
  meanFun <- function(covarVec, ranEffect) {
    sum(fixedCoefs*c(1,covarVec)) + ranEffect
  }
  covValuesWithIntercept <- cbind(Intercept = 1, covValues)
  lapplyInternal <- function(obsIndex) {
    derivMean <- c(1, x, rep(1, length(y)))
    (y[obsIndex] - meanFun(covValuesWithIntercept[obsIndex, ], ranEffects[obsIndex])) * derivMean
  }
  allDerivs <- lapply(seq_along(y), lapplyInternal)
  -do.call('+', allDerivs)/MEvar + c(fixedCoefs, ranEffects) %*% MRAprecision
}

.FisherInfo <- function(covValues, MEvar, MRAprecision) {
  derivVecs <- lapply(1:nrow(covValues), function(obsIndex) c(1, covValues[obsIndex, ], rep(1, nrow(covValues))))
  derivMat <- do.call("cbind", derivVecs)
  MEvar * t(derivMat) %*% derivMat + MRAprecision
}

.FisherScoringAlgo <- function(tol = 1e-5, startingVal = rep(0, nrow(invFisherInfo)), invFisherInfo, scoreFun) {
  updatedVal <- startingVal + invFisherInfo %*% scoreFun(startingVal)
  if (max(updatedVal - startingVal) > tol) {
    updatedVal <- .FisherScoringAlgo(tol = tol, startingVal = updatedVal, invFisherInfo = invFisherInfo, scoreFun = scoreFun)
  }
  updatedVal
}

# In our situation, the curvature is constant, hence the absence of a position argument.

.curvatureFun <- function(covValues, MEvar, MRAprecision) {
  curvatures <- lapply(1:nrow(covValues), function(rowIndex) c(1, covValues[rowIndex, ], rep(1, nrow(covValues)))^2)
  do.call("+", curvatures)/MEvar + t(rowSums(MRAprecision))
}

.logGMRFprior <- function(paraValuesList, SigmaValuesList) {
  sum(mapply(paraValues = paraValuesList, varPar = SigmaValuesList, FUN = function(paraValues, varPar) {
    if (!is.matrix(varPar)) {
      varPar <- varPar*diag(length(paraValues))
    }
    mvtnorm::dmvnorm(x = paraValues, sigma = varPar, log = TRUE)
  }))
}

computeLogLik <- function(gridObj, covFct, fixedEffectParVec = NULL) {
  .computeWmats(gridObj = gridObj, covFct = covFct)
  .setBtips(gridObj)
  .setSigmaTips(gridObj)
  .setAtildeTips(gridObj)
  .recurseA(gridObj)
  .setOmegaTildeTips(gridObj, fixedEffectVec = fixedEffectParVec)
  .recurseOmega(gridObj)

  .computeUtips(gridObj, fixedEffectVec = fixedEffectParVec)
  .recurseU(gridObj)

  .computeDtips(gridObj)
  .recurseD(gridObj)

  gridObj$logLik <- -(gridObj$d + gridObj$u)/2
  invisible()
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
    x$SigmaInverse <- Matrix::chol2inv(Matrix::chol(x$Sigma))
    invisible()
  }) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.setAtildeTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  M <- gridObj$M

  modifyTip <- function(x) {

    x$Atilde <- vector('list', x$depth+1) # First level of the hierarchy is for k in Atilde^{k,l}.
    x$Atilde <- lapply(0:x$depth, function(k) vector('list', x$depth)) # Second level is for l in Atilde^{k,l}, same order as for k.
    indexGrid <- expand.grid(k = 0:(M-1), l = 0:(M-1)) # We don't need A^{M,M}_{j_1, ..., j_M}
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
    lapply(brickObj$childBricks, .recurseD)
  }

  brickObj$d <- log(det(brickObj$KtildeInverse)) - log(det(brickObj$Kinverse)) + sum(sapply(brickObj$childBricks, FUN = function(childBrick) childBrick$d))
  invisible()
}

# If elements of fixedEffectVec and the columns in "dataset" element are named, the function will try to match them. The name of the intercept term does not matter.
# If the names of the covariates don't match, an error is produced. If names are missing in either, it is assumed that the first term in fixedEffectVec is for the intercept, and that the order of the other terms matches that of the columns in dataset@data (the response can be in any column, as long as the ordering of the covariate columns is the same as in fixedEffectVec).
# If the response column in dataset@data is not called "y", it is assumed that the first column is the response.

.setOmegaTildeTips <- function(gridObj, fixedEffectVec) {
  allTips <- .tipAddresses(gridObj)

  modifyTip <- function(x) {
    x$omegaTilde <- replicate(n = x$depth+1, expr = vector('list', x$depth + 1), simplify = FALSE)

    rescaledObservations <- .rescaleObservations(gridObj$dataset@data[x$observations, ], fixedEffectVec)

    for (k in 0:(gridObj$M-1)) {
      x$omegaTilde[[k+1]] <- t(x$BmatList[[k+1]]) %*% x$SigmaInverse %*% rescaledObservations
    }
    invisible()
  }
  lapply(allTips, FUN = modifyTip) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.rescaleObservations <- function(obsData, fixedEffectVec) {
  responseCol <- match("y", colnames(obsData)) # Will return NA if obsData does not have colnames.
  if (is.na(responseCol)) {
    responseCol <- 1
  }
  covData <- as.matrix(cbind(Placeholderz = 1, obsData[ , -responseCol, drop = FALSE])) # The additional column is for the intercept.

  if (!is.null(names(fixedEffectVec)) & !is.null(colnames(obsData))) {
    matchingOrder <- match(names(fixedEffectVec), colnames(covData))
    matchingOrder <- c(1, matchingOrder[!is.na(matchingOrder)]) # The NA appears because of the intercept term.
    covData <- covData[ , matchingOrder]
  }
  rescaledObservations <- obsData[ , responseCol]
  if  (!is.null(fixedEffectVec)) {
    rescaledObservations <- rescaledObservations - drop(covData %*% fixedEffectVec)
  }
  unname(rescaledObservations)
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
  })

 brickObj$omega <- omegaMatrices

  brickObj$omegaTilde <- lapply(0:m, FUN = function(k) {
    brickObj$omega[[k+1]] - brickObj$A[[k+1]][[m+1]] %*% brickObj$Ktilde %*% brickObj$omega[[m+1]]
  })
  invisible()
}

.computeUtips <- function(gridObj, fixedEffectVec) {
  allTips <- .tipAddresses(gridObj)

  modifyTip <- function(x) {
    rescaledObservations <- .rescaleObservations(gridObj$dataset@data[x$observations, ], fixedEffectVec)
    x$u <- t(rescaledObservations) %*% x$SigmaInverse %*% rescaledObservations ## QUADRATIC FORM: MAY BE OPTIMIZED.
    invisible()
  }
  lapply(allTips, FUN = modifyTip)
  invisible()
}

.recurseU <- function(brickObj) {
  if (is.null(brickObj$childBricks)) {
    return(invisible())
  }

  if (is.null(brickObj$childBricks[[1]]$u)) {
    lapply(brickObj$childBricks, .recurseU)
  }

  brickObj$u <- -t(brickObj$omega[[brickObj$depth+1]]) %*% brickObj$Ktilde %*% brickObj$omega[[brickObj$depth + 1]] + sum(sapply(brickObj$childBricks, FUN = function(childBrick) childBrick$u)) ## QUADRATIC FORM, MIGHT BE POSSIBLE TO OPTIMIZE.
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

# Prediction functions

predict.Spacetimegrid <- function(gridObj, spacetimeCoor) {
  .computeBtildeTips(gridObj, spacetimeCoor)
  .computeLtips(gridObj, spacetimeCoor)
  .computeVtips(gridObj, spacetimeCoor)
  meanVector <- .computeMeanDist(gridObj, spacetimeCoor)
  covMatrix <- .computeCovMat(gridObj, spacetimeCoor)
  list(predictions = meanVector, covMat = covMatrix)
}

.computeBtildeTips <- function(gridObj, locations) {
  allTips <- .tipAddresses(gridObj)
  # We'll fit Btilde^{l,k} in a list of lists

  # First for l = M... We need the tail elements only. We need b_{j_1, ..., j_{M-1}}(predict locations).
  lapply(allTips, .computeBtildeTipsInternal, locations = locations)
}

.computeBtildeInternal <- function(brickObj, locations) {
  BtildeList <- vector(mode = 'list', length = gridObj$M)
  BtildeList <- lapply(seq_along(BtildeList), function(index) vector(mode = 'list', length = index))
  gridObj <- .getTopEnvirAddress(brickObj)
  vZero <- gridObj$covFct(locations, gridObj$knotPositions)

  updatedV <- vZero
  pointedBrick <- gridObj

  for (i in 1:(gridObj$M-1)) {
    lapply(pointedBrick$childBricks, function(brickPointer) {

    })
    updatedV <- updatedV -
  }
}
