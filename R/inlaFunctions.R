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
  gridObj$covFct <- covFct
  .computeWmats(gridObj)
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

  gridObj$logLik <- drop(-(gridObj$d + gridObj$u)/2)
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

.computeWmats <- function(gridObj) {
  gridObj$Kinverse <- gridObj$covFct(gridObj$knotPositions, gridObj$knotPositions)
  gridObj$K <- Matrix::chol2inv(Matrix::chol(gridObj$Kinverse))
  gridObj$WmatList <- list(gridObj$K)

  .computeWchildBrick <- function(brickObj) {
    brickObj$WmatList <- vector("list", brickObj$depth+1) # Depth is from 0 to M.
    brickList <- .getAllParentAddresses(brickObj, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    brickObj$WmatList[[1]] <- gridObj$covFct(brickObj$knotPositions, brickList[[1]]$knotPositions)
    m <- brickObj$depth

    for (l in 1:m) { # First element in brickList is the top brick, i.e. at resolution 0.
      firstTerm <- gridObj$covFct(brickObj$knotPositions, brickList[[l+1]]$knotPositions)
      makeSumTerm <- function(k) {
        brickObj$WmatList[[k+1]] %*% brickList[[k+1]]$K %*% t(brickList[[l+1]]$WmatList[[k+1]])
      }
      secondTerm <- Reduce("+", lapply(0:(l-1), makeSumTerm))

      brickObj$WmatList[[l+1]] <- firstTerm - secondTerm
    }

    if (brickObj$depth < gridObj$M) {
      brickObj$Kinverse <- tail(brickObj$WmatList, n = 1)[[1]]
      brickObj$K <- Matrix::chol2inv(Matrix::chol(brickObj$Kinverse))
    }

    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, .computeWchildBrick, covFct = covFct)
    }
    invisible()
  }

  if (!is.null(gridObj$childBricks)) {
    lapply(gridObj$childBricks, .computeWchildBrick, covFct = gridObj$covFct)
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
  if (is.null(gridObj$logLik)) {
    stop("Perform likelihood evaluation first. \n")
  }
  # .computeLmatrices(gridObj, spacetimeCoor)
  .computeUpredictTips(gridObj, spacetimeCoor)
  .computeVpredictTips(gridObj, spacetimeCoor)
  .computeBtildeTips(gridObj, spacetimeCoor)

  meanVector <- .computeMeanValues(gridObj, spacetimeCoor)
  varianceEst <- .computeVar(gridObj, spacetimeCoor)
  list(predictions = meanVector, variance = covMatrix)
}

.computeBtildeTips <- function(gridObj, locations) {
  .computeBknots(gridObj)
  .computeBpred(gridObj, locations)

  allTips <- .tipAddresses(gridObj)
  .computeBtildeInternal <- function(tipAddress) {
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    lapply(seq_along(tipAddress$Btilde), function(index) tipAddress$Btilde <- vector('list', length = index)) # The index k can take values 0 to l-1.

    funForDefinition <- function(l,k) {
      tipAddress$Btilde[[l+1]][[k+1]] <- tipAddress$bPred[[k+1]] - tipAddress$UpredList[[gridObj$M+1]] %*% brickList[[l+1]]$SigmaInverse %*% brickList[[l+1]]$WmatList[[k+1]]
      invisible()
    }
    lapply((gridObj$M-1):(gridObj$M-3), funForDefinition, l = gridObj$M)
    # We can now launch the recursion.
    funForRecursion <- function(l,k) {
      tipAddress$Btilde[[l+1]][[k+1]] <- tipAddress$Btilde[[l+2]][[k+1]] - tipAddress$Btilde[[l+2]][[l+1]] %*% brickList[[l+1]]$Ktilde %*% brickList[[l+1]]$A[[l+1]][[k+1]]
      invisible()
    }
    for (m in (gridObj$M-1):1) {
      funForRecursion(m, m-2)
      funForRecursion(m, m-1)
      funForRecursion(m-1, m-2)
    }
    invisible()
  }
  lapply(allTips, .computeBtildeInternal)
}

.computeBknots <- function(gridObj) {
  initiateBs <- function(brickObj) {
    brickObj$bKnots <- vector('list', length = brickObj$depth)
    brickObj$bKnots[[1]] <- gridObj$covFct(brickObj$knotPositions, gridObj$knotPositions)
    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, childFun)
    }
    invisible()
  }
  lapply(gridObj$childBricks, initiateBs)

  completeBknotsLevel <- function(brickObj) {
    firstUnresolved <- match(TRUE, sapply(brickObj$bKnots, function(x) is.null(x)))
    if (is.na(firstUnresolved)) {
      return(invisible())
    }
    knots1 <- brickObj$knotPositions
    brickList <- .getAllParentAddresses(brickObj = brickObj, flip = TRUE)
    knots2 <- parent.env(brickObj)$knotPositions
    currentV <- gridObj$covFct(knots1, knots2)
    for (i in 2:firstUnresolved) {
      currentV <- currentV - t(brickObj$bKnots[[i-1]]) %*% brickList[[i-1]]$K %*% parent.env(brickObj)$bKnots[[i-1]]
    }
    brickObj$bKnots[[firstUnresolved]] <- currentV
    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, completeBknots)
    }
    invisible()
  }
  lapply(1:gridObj$M, function(levelIndex) {
    lapply(gridObj$childBricks, completeBknotsLevel)
  })
  invisible()
}

.computeBpreds <- function(gridObj, locations) {

  allTips <- .tipAddresses(gridObj)

  recurseBpred <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    tipAddress$bPred <- vector('list', length = gridObj$M)
    tipAddress$bPred[[1]] <- gridObj$covFct(subLocations, gridObj$knotPositions)

    brickList <- .getAllParentAddresses(brickObj = brickObj, flip = TRUE)

    recurseV <- function(locations1, locations2, level) {

      currentV <- gridObj$covFct(locations1, locations2) - t(tipAddress$bPred[[1]]) %*% gridObj$K %*% brickList[[level]]$bKnots[[1]]
      if (level == 1) {
        tipAddress$bPred[[2]] <- currentV
        return(invisible())
      }

      for (i in 1:(level-1)) {
        currentV <- currentV - t(tipAddress$bPred[[i+1]]) %*% brickList[[i+1]]$K %*% brickList[[level]]$bKnots[[i+1]]
      }
      currentV
    }

    for (i in 2:(gridObj$M-1)) {
      tipAddress$bPred[[i]] <- recurseV(subLocations, brickList[[i]], level = i-1) - t(tipAddress$bPred[[i-1]]) %*% brickList[[i-1]]$K %*% brickList[[i]]$bKnots[[i-1]]
    }
    invisible()
  }
  lapply(allTips, recurseBpred)
}

# .computeLmatrices <- function(gridObj, locations) {
#   allTips <- .tipAddresses(gridObj)
#   computeTipLevel <- function(tipAddress) {
#     subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
#     brickList <- .getAllParentAddresses(brickObj = brickObj, flip = TRUE)
#     currentV <- gridObj$covFct(subLocations, gridObj$knotPositions) - tipAddress$bPred[[1]] %*% gridObj$K %*% tipAddress$WmatList[[1]]
#     if (gridObj$M < 2) {
#       return(currentV)
#     }
#     for (i in 2:gridObj$M) {
#       currentV <- currentV - tipAddress$bPred[[i]] %*% brickList[[i]]$K %*% tipAddress$WmatList[[i]]
#     }
#     tipAddress$L <- currentV
#     invisible()
#   }
#   lapply(allTips, computeTipLevel)
# }

.computeMeanValues <- function(gridObj, locations) {
  allTips <- .tipAddresses(gridObj)

  predictMeanInTip <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    computeSumTerm <- function(m) {
      tipAddress$Btilde[[m+2]][[m+1]] %*% brickList[[m+1]]$Ktilde %*% brickList[[m+1]]$omega
    }
    firstTerm <- do.call('+', lapply(0:gridObj$M-1, computeSumTerm))
    secondTerm <- tipAddress$UpredList[[gridObj$M+1]] %*% tipAddress$SigmaInverse %*% gridObj@dataset[tipAddress$observations, ]
    predictionMean <- firstTerm + secondTerm
    names(predictionMean) <- rownames(subLocations@data)
    predictionMean
  }
  unorderedPredictions <- lapply(allTips, predictMeanInTip)
  do.call('c', unorderedPredictions)[rownames(locations@data)]
}

.computeVar <- function(gridObj, locations) {
  allTips <- .tipAddresses(gridObj)
  variancesInTip <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    computeSumTermVar <- function(m) {
      diag(tipAddress$Btilde[[m+2]][[m+1]])^2 * diag(brickList[[m+1]]$Ktilde)
    }
    firstTerm <- do.call('+', lapply(0:(gridObj$M-1), computeSumTermVar))
    secondTerm <- tipAddress$VpredMat - tipAddress$UpredList[[gridObj$M+1]] %*% tipAddress$SigmaInverse %*% t(tipAddress$UpredList[[gridObj$M+1]])
    variances <- firstTerm + secondTerm
    names(variances) <- rownames(subLocations@data)
    variances
  }
  predictVars <- do.call('c', lapply(allTips, variancesInTip))
  predictVars[rownames(locations@data)]
}

.computeUpredictTips <- function(gridObj, locations) {
  allTips <- .tipAddresses(gridObj)
  computeUlistTips <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    tipAddress$UpredList <- vector('list', length = gridObj$M+1)
    tipAddress$UpredList[[1]] <- gridObj$covFct(subLocations, gridObj$knotPositions)

    for (l in 1:gridObj$M) {
      firstTerm <- gridObj$covFct(subLocations, brickList[[l+1]]$knotPositions)
      secondTermList <- lapply(0:(l-1), function(k) tipAddress$UpredList[[k+1]] %*% brickList[[k+1]]$K %*% brickList[[l+1]]$W[[k+1]])

      tipAddress$UpredList[[l+1]] <- firstTerm - Reduce(f = '+', x = secondTermList)
    }
    invisible()
  }
  lapply(allTips, computeUlistTips)
  invisible()
}

.computeVpredictTips <- function(gridObj, locations) {
  allTips <- .tipAddresses(gridObj)
  computeVpredictMat <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    firstTerm <- gridObj$covFct(subLocations, subLocations)
    secondTermList <- lapply(0:(gridObj$M-1), function(k) tipAddress$UpredList[[k+1]] %*% brickList[[k+1]]$K %*% t(tipAddress$UpredList[[k+1]]))
    tipAddress$VpredMat <- firstTerm - Reduce('+', secondTermList)
    invisible()
  }
  lapply(allTips, computeVpredictMat)
  invisible()
}
