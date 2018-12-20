computeLogLik <- function(gridObj, covFct, fixedEffectParVec = NULL, cl = NULL) {
  gridObj$covFct <- covFct
  .computeWmats(gridObj, cl)
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
  .clean(gridObj)
  invisible()
}

.computeWmats <- function(gridObj, cl = NULL) {
  gridObj$Kinverse <- gridObj$covFct(gridObj$knotPositions, gridObj$knotPositions)
  # gridObj$K <- Matrix::chol2inv(Matrix::chol(gridObj$Kinverse))
  gridObj$K <- solve(gridObj$Kinverse)
  gridObj$WmatList <- list(gridObj$Kinverse)

  computeWchildBrick <- function(brickObj) {
    if (is.null(brickObj$knotPositions)) {
      brickObj$WmatList <- NULL
      brickObj$Kinverse <- NULL
      brickObj$K <- NULL
      return(invisible())
    }
    brickObj$WmatList <- vector("list", brickObj$depth+1) # Depth is from 0 to M.
    brickList <- .getAllParentAddresses(brickObj, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    brickObj$WmatList[[1]] <- gridObj$covFct(brickObj$knotPositions, gridObj$knotPositions)
    m <- brickObj$depth

    for (l in 1:m) {
      firstTerm <- gridObj$covFct(brickObj$knotPositions, brickList[[l+1]]$knotPositions)
      makeSumTerm <- function(k) {
        brickObj$WmatList[[k+1]] %*% brickList[[k+1]]$K %*% t(brickList[[l+1]]$WmatList[[k+1]])
      }
      secondTerm <- Reduce("+", lapply(0:(l-1), makeSumTerm))

      brickObj$WmatList[[l+1]] <- firstTerm - secondTerm
    }

    if (m < gridObj$M) {
      brickObj$Kinverse <-brickObj$WmatList[[m+1]]
      # brickObj$K <- Matrix::chol2inv(Matrix::chol(brickObj$Kinverse))
      brickObj$K <- solve(brickObj$Kinverse)
    }

    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, computeWchildBrick)
    }
    brickObj
  }

  if (!is.null(gridObj$childBricks)) {
    if (is.null(cl)) {
      lapply(gridObj$childBricks, computeWchildBrick)
    } else {
      newBricks <- parallel::parLapply(cl = cl, X = gridObj$childBricks, fun = computeWchildBrick)
      gridObj$childBricks <- newBricks
    }
  }
  invisible()
}

.setBtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  lapply(allTips, FUN = function(x) {
    x$BmatList <- x$WmatList[1:gridObj$M]
    invisible()
  })
  invisible()
}

.setSigmaTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  lapply(allTips, FUN = function(x) {
    if (!is.null(x$WmatList)) {
      x$Sigma <- x$WmatList[[gridObj$M+1]]
      # x$SigmaInverse <- Matrix::chol2inv(Matrix::chol(x$Sigma))
      x$SigmaInverse <- solve(x$Sigma)
    }
    invisible()
  }) # The last element in Bmat will be NULL, since B^M_{j_1, ..., j_M} is undefined,
  invisible()
}

.setAtildeTips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)
  M <- gridObj$M

  indexGrid <- expand.grid(k = 0:(M-1), l = 0:(M-1)) # We don't need A^{M,M}_{j_1, ..., j_M}.
  indexGrid <- subset(indexGrid, subset = k >= l)

  modifyTip <- function(tipAddress) {
    if (is.null(tipAddress$knotPositions)) {
      tipAddress$Atilde <- NULL
      return(invisible())
    }
    tipAddress$Atilde <- vector('list', tipAddress$depth+1) # First level of the hierarchy is for k in Atilde^{k,l}.
    tipAddress$Atilde <- lapply(0:tipAddress$depth, function(k) {
      tipAddress$Atilde[[k+1]] <- vector('list', length = k+1) # Second level is for l in Atilde^{k,l}.
    })

    for (i in 1:nrow(indexGrid)) {
      k <- indexGrid[i, 1]
      l <- indexGrid[i, 2]
      tipAddress$Atilde[[k + 1]][[l + 1]] <- t(tipAddress$BmatList[[k + 1]]) %*% tipAddress$SigmaInverse %*% tipAddress$BmatList[[l + 1]]
      # tipAddress$Atilde[[l + 1]][[k + 1]] <- t(tipAddress$Atilde[[k + 1]][[l + 1]])
    }
    invisible()
  }
  lapply(allTips, FUN = modifyTip)
  invisible()
}

.recurseA <- function(brickObj) {
  if (is.null(brickObj$childBricks)) {
    return(NULL)
  }

  if (is.null(brickObj$childBricks[[1]]$Atilde)) {
    lapply(brickObj$childBricks, .recurseA)
  }

  if (is.null(brickObj$knotPositions)) {
    brickObj$Atilde <- NULL
    brickObj$Ktilde <- NULL
    brickObj$KtildeInverse <- NULL
    brickObj$A <- NULL
    return(invisible())
  }

  m <- brickObj$depth
  indexGrid <- expand.grid(k = 0:m, l = 0:m)
  indexGrid <- subset(indexGrid, subset = k >= l)

  mapplyFunction <- function(k,l) { # +1 to adjust indices for R.
    Reduce("+", Filter(f = function(x) !is.null(x), x = lapply(brickObj$childBricks, function(x) {
      x$Atilde[[k + 1]][[l + 1]]
    })))
  }
  Amatrices <- mapply(k = indexGrid$k, l = indexGrid$l, FUN = mapplyFunction, SIMPLIFY = FALSE)
  brickObj$A <- .reformatList(matrixList = Amatrices, indexGrid = indexGrid)

  brickObj$KtildeInverse <- brickObj$Kinverse + brickObj$A[[m + 1]][[m + 1]] # +1 to adjust indices for R.
  # brickObj$Ktilde <- Matrix::chol2inv(chol(brickObj$KtildeInverse))
  brickObj$Ktilde <- solve(brickObj$KtildeInverse)
  mapplyFunction <- function(k, l) {
    brickObj$A[[k + 1]][[l + 1]] - brickObj$A[[k + 1]][[m + 1]] %*% brickObj$Ktilde %*% brickObj$A[[m + 1]][[l + 1]] # m+1 to adjust indices for R.
  }

  AtildeMatrices <- mapply(k = indexGrid$k, l = indexGrid$l, FUN = mapplyFunction, SIMPLIFY = FALSE)
  brickObj$Atilde <- .reformatList(matrixList = AtildeMatrices, indexGrid = indexGrid)

  invisible()
}

.computeDtips <- function(gridObj) {
  allTips <- .tipAddresses(gridObj)

  modifyTip <- function(tipAddress) {
    if (is.null(tipAddress$Sigma)) {
      tipAddress$d <- NULL
      return(invisible())
    }
    tipAddress$d <- determinant(tipAddress$Sigma, logarithm = TRUE)$modulus
    attributes(tipAddress$d) <- NULL # I'm stripping the attributes of the object returned by determinant.
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
  firstTerm <- determinant(brickObj$KtildeInverse, logarithm = TRUE)$modulus
  secondTerm <- determinant(brickObj$Kinverse, logarithm = TRUE)$modulus
  attributes(firstTerm) <- attributes(secondTerm) <- NULL
  brickObj$d <- firstTerm - secondTerm + sum(do.call('c', Filter(f = function(x) !is.null(x), x = lapply(brickObj$childBricks, FUN = function(childBrick) childBrick$d))))
  if (brickObj$d == Inf) browser()
  invisible()
}

# If elements of fixedEffectVec and the columns in "dataset" element are named, the function will try to match them. The name of the intercept term does not matter.
# If the names of the covariates don't match, an error is produced. If names are missing in either, it is assumed that the first term in fixedEffectVec is for the intercept, and that the order of the other terms matches that of the columns in dataset@data (the response can be in any column, as long as the ordering of the covariate columns is the same as in fixedEffectVec).
# If the response column in dataset@data is not called "y", it is assumed that the first column is the response.

.setOmegaTildeTips <- function(gridObj, fixedEffectVec) {
  allTips <- .tipAddresses(gridObj)

  modifyTip <- function(tipAddress) {
    if (is.null(tipAddress$knotPositions)) {
      tipAddress$omegaTilde <- NULL
      return(invisible())
    }
    tipAddress$omegaTilde <- replicate(n = tipAddress$depth + 1, expr = vector('list', tipAddress$depth + 1), simplify = FALSE) ## Is this correct? omegaTilde will be a list of lists although it has only one index. The following loop might even overwrite this step.

    # rescaledObservations <- .rescaleObservations(gridObj$dataset@data[tipAddress$observations, ], fixedEffectVec)

    for (k in 0:(gridObj$M - 1)) {
      tipAddress$omegaTilde[[k + 1]] <- t(tipAddress$BmatList[[k + 1]]) %*% tipAddress$SigmaInverse %*% gridObj$dataset[tipAddress$observations]@data[ , "y"]
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
    Reduce('+', Filter(f = function(x) !is.null(x), x = lapply(brickObj$childBricks, function(x) {
      x$omegaTilde[[k+1]]
    })))
  })

  brickObj$omega <- omegaMatrices

  brickObj$omegaTilde <- lapply(0:m, FUN = function(k) {
    brickObj$omega[[k+1]] - brickObj$A[[k+1]][[m+1]] %*% brickObj$Ktilde %*% brickObj$omega[[m+1]]
  })
  invisible()
}

.computeUtips <- function(gridObj, fixedEffectVec) {
  allTips <- .tipAddresses(gridObj)

  modifyTip <- function(tipAddress) {
    # rescaledObservations <- .rescaleObservations(gridObj$dataset@data[tipAddress$observations, ], fixedEffectVec)
    if (is.null(tipAddress$knotPositions)) {
      tipAddress$u <- NULL
      return(invisible())
    }
    tipAddress$u <- drop(t(gridObj$dataset[tipAddress$observations]@data[, "y"]) %*% tipAddress$SigmaInverse %*% gridObj$dataset[tipAddress$observations]@data[, "y"]) ## QUADRATIC FORM: MAY BE OPTIMIZED.
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

  brickObj$u <- -t(brickObj$omega[[brickObj$depth+1]]) %*% brickObj$Ktilde %*% brickObj$omega[[brickObj$depth + 1]] + sum(do.call('c', Filter(f = function(x) !is.null(x), x = lapply(brickObj$childBricks, FUN = function(childBrick) childBrick$u)))) ## QUADRATIC FORM, MIGHT BE POSSIBLE TO OPTIMIZE.
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

predict.Spacetimegrid <- function(gridObj, spacetimeCoor, cl = NULL) {
  if (is.null(gridObj$logLik)) {
    stop("Perform likelihood evaluation first. \n")
  }
  # .computeLmatrices(gridObj, spacetimeCoor)
  .computeUpredictTips(gridObj, spacetimeCoor, cl)
  .computeVpredictTips(gridObj, spacetimeCoor, cl)
  .computeBtildeTips(gridObj, spacetimeCoor, cl)

  meanVector <- .computeMeanValues(gridObj, spacetimeCoor)
  varianceEst <- .computeVar(gridObj, spacetimeCoor)
  list(predictions = meanVector, variance = varianceEst)
}

.computeBtildeTips <- function(gridObj, locations, cl) {
  .computeBknots(gridObj)
  .computeBpreds(gridObj, locations)

  allTips <- .tipAddresses(gridObj)
  .computeBtildeInternal <- function(tipAddress) {
    if (is.null(tipAddress$bPred)) {
      tipAddress$Btilde <- NULL
      return(invisible())
    }
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    lapply(seq_along(tipAddress$Btilde), function(index) tipAddress$Btilde <- vector('list', length = index)) # The index k can take values 0 to l-1.

    funForDefinition <- function(l,k) {
      tipAddress$Btilde[[l+1]][[k+1]] <- tipAddress$bPred[[k+1]] - tipAddress$UpredList[[gridObj$M+1]] %*% brickList[[l+1]]$SigmaInverse %*% brickList[[l+1]]$WmatList[[k+1]]
      invisible()
    }
    lapply((gridObj$M-1):max(gridObj$M-3, 0), funForDefinition, l = gridObj$M)

    # We can now launch the recursion.

    funForRecursion <- function(l,k) {
      tipAddress$Btilde[[l+1]][[k+1]] <- tipAddress$Btilde[[l+2]][[k+1]] - tipAddress$Btilde[[l+2]][[l+1]] %*% brickList[[l+1]]$Ktilde %*% brickList[[l+1]]$A[[l+1]][[k+1]]
      invisible()
    }
    for (k in (gridObj$M-1):1) { # What happens when m is 1?
      funForRecursion(k, k-1)
      if (k > 1) funForRecursion(k, k-2)
      if (k > 2) funForRecursion(k, k-3)
    }
    invisible()
  }
  if (is.null(cl)) {
    lapply(allTips, .computeBtildeInternal)
  } else {
    parallel::parLapply(cl = cl, X = allTips, fun = .computeBtildeInternal)
  }
  invisible()
}

.computeBknots <- function(gridObj) {
  initiateBs <- function(brickObj) {
    brickObj$bKnots <- vector('list', length = brickObj$depth)
    brickObj$bKnots[[1]] <- gridObj$covFct(brickObj$knotPositions, gridObj$knotPositions)
    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, initiateBs)
    }
    invisible()
  }
  lapply(gridObj$childBricks, initiateBs)

  completeBknotsLevel <- function(brickObj, level) {
    knots1 <- brickObj$knotPositions
    brickList <- .getAllParentAddresses(brickObj = brickObj, flip = FALSE)
    knots2 <- brickList[[level+1]]$knotPositions
    currentV <- gridObj$covFct(knots1, knots2)

    for (i in 0:(level-1)) {
      currentV <- currentV - brickObj$bKnots[[i+1]] %*% brickList[[i+1]]$K %*% t(brickList[[level+1]]$bKnots[[i+1]]) ### MOVED THE t(...) TO THE RIGHT SIDE OF K. This differs from the paper, but solves the problem with the matrix dimensions.
    }
    brickObj$bKnots[[level+1]] <- currentV
    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, completeBknotsLevel, level = level)
    }
    invisible()
  }
  lapply(1:(gridObj$M-1), function(levelIndex) {
    lapply(getLayer(gridObj, m = levelIndex+1), completeBknotsLevel, level = levelIndex)
  })
  invisible()
}

.computeBpreds <- function(gridObj, locations) {

  allTips <- .tipAddresses(gridObj)

  recurseBpred <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    tipAddress$bPred <- vector('list', length = gridObj$M+1)
    if (length(subLocations) == 0) {
      tipAddress$bPred <- NULL
      return(invisible())
    }
    tipAddress$bPred[[1]] <- gridObj$covFct(subLocations, gridObj$knotPositions)

    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) ## Lists bricks from root to current tip (included).

    recurseBpred <- function(level) {
      currentV <- gridObj$covFct(subLocations, brickList[[level+1]]$knotPositions) - tipAddress$bPred[[1]] %*% gridObj$K %*% t(brickList[[level+1]]$bKnots[[1]]) ### SAME ISSUE AS IN .computeBknots with the position of the t(...). ###
      if (level == 1) {
        tipAddress$bPred[[2]] <- currentV
        return(invisible())
      }

      for (i in 2:level) {
        currentV <- currentV - tipAddress$bPred[[i]] %*% brickList[[i]]$K %*% t(brickList[[level+1]]$bKnots[[i]]) # Changed the position of the t to ensure correspondence.
      }
      tipAddress$bPred[[level+1]] <- currentV
      invisible()
    }

    for (levelIndex in 1:gridObj$M) {
      recurseBpred(levelIndex)
    }
    invisible()
  }
  lapply(allTips, recurseBpred)
}

.computeMeanValues <- function(gridObj, locations) {
  allTips <- .tipAddresses(gridObj)

  predictMeanInTip <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    if (length(subLocations) == 0) {
      return(NULL)
    }
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    computeSumTerm <- function(m) {
      drop(tipAddress$Btilde[[m+2]][[m+1]] %*% brickList[[m+1]]$Ktilde %*% brickList[[m+1]]$omega[[m+1]])
    }
    firstTerm <- Reduce('+', lapply(0:(gridObj$M-1), computeSumTerm))
    secondTerm <- drop(tipAddress$UpredList[[gridObj$M+1]] %*% tipAddress$SigmaInverse %*% gridObj$dataset[tipAddress$observations]@data[ , 1])
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
    if (length(subLocations) == 0) {
      return(NULL)
    }
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    computeSumTermCov <- function(m) {
      tipAddress$Btilde[[m+2]][[m+1]] %*% brickList[[m+1]]$Ktilde %*% t(tipAddress$Btilde[[m+2]][[m+1]])
    }
    firstTerm <- Reduce('+', lapply(0:(gridObj$M-1), computeSumTermCov))
    secondTerm <- tipAddress$VpredMat - tipAddress$UpredList[[gridObj$M+1]] %*% tipAddress$SigmaInverse %*% t(tipAddress$UpredList[[gridObj$M+1]])
    covariances <- firstTerm + secondTerm
    rownames(covariances) <- colnames(covariances) <- rownames(subLocations@data)
    diag(covariances) # Covariances could be returned too, but the code for now only returns variances.
  }
  predictVars <- do.call('c', lapply(allTips, variancesInTip))
  predictVars[rownames(locations@data)]
}

.computeUpredictTips <- function(gridObj, locations, cl = NULL) {
  allTips <- .tipAddresses(gridObj)
  computeUlistTips <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    if (length(subLocations) == 0) {
      tipAddress$UpredList <- NULL
      return(invisible())
    }
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    tipAddress$UpredList <- vector('list', length = gridObj$M+1)
    tipAddress$UpredList[[1]] <- gridObj$covFct(subLocations, gridObj$knotPositions)

    for (l in 1:gridObj$M) {
      firstTerm <- gridObj$covFct(subLocations, brickList[[l+1]]$knotPositions)
      secondTermList <- lapply(0:(l-1), function(k) tipAddress$UpredList[[k+1]] %*% brickList[[k+1]]$K %*% t(brickList[[l+1]]$WmatList[[k+1]]))

      tipAddress$UpredList[[l+1]] <- firstTerm - Reduce(f = '+', x = secondTermList)
    }
    invisible()
  }
  if (is.null(cl)) {
    lapply(allTips, computeUlistTips)
  } else {
    parallel::parLapply(cl = cl, X = allTips, fun = computeUlistTips)
  }
  invisible()
}

.computeVpredictTips <- function(gridObj, locations, cl = NULL) {
  allTips <- .tipAddresses(gridObj)
  computeVpredictMat <- function(tipAddress) {
    subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
    if (length(subLocations) == 0) {
      tipAddress$VpredMat <- NULL
      return(invisible())
    }
    brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) # Lists parent nodes in order from root to itself (included).
    firstTerm <- gridObj$covFct(subLocations, subLocations)
    secondTermList <- lapply(0:(gridObj$M-1), function(k) tipAddress$UpredList[[k+1]] %*% brickList[[k+1]]$K %*% t(tipAddress$UpredList[[k+1]]))
    tipAddress$VpredMat <- firstTerm - Reduce('+', secondTermList)
    invisible()
  }
  if (is.null(cl)) {
    lapply(allTips, computeVpredictMat)
  } else {
    parallel::parLapply(cl = cl, X = allTips, fun = computeVpredictMat)
  }
  invisible()
}

.clean <- function(gridObj) {
  rm(list = setdiff(ls(gridObj), c("knotPositions", "observations", "covFct", "dimensions", "logLik", "depth", "M", "childBricks", "dataset", "WmatList", "K", "Ktilde", "Sigma", "SigmaInverse", "A", "omega")), envir = gridObj)
  internalFun <- function(brickObj) {
    rm(list = setdiff(ls(brickObj), c("knotPositions", "observations", "dimensions", "Ktilde", "depth", "childBricks", "WmatList", "K", "Sigma", "SigmaInverse", "A", "omega")), envir = brickObj)
    if (!is.null(brickObj$childBricks)) {
      lapply(brickObj$childBricks, internalFun)
    }
    invisible()
  }
  lapply(gridObj$childBricks, internalFun)
  invisible()
}

predictMRA <- function(treePointer, predictLocations, timePoints, covariatesForPreds) {
  predictions <- predictMRArcpp(treePointer = treePointer, predSpatialCoor = predictLocations, predTime = timePoints, covariatesAtPredLocs = covariatesForPreds)
  predictions$means <- do.call("c", predictions$means)
  predictions$covar <- lapply(predictions$covMatricesNoDim, function(x) {
    numRows <- numCols <- sqrt(length(x))
    dim(x) <- c(numRows, numCols)
    x
  })
  predictions <- predictions[-which(names(predictions) == "covMatricesNoDim")]
  predictions
}

bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if (!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if (N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M - 1L), nrow = k)[, rep(seq_len(N), each = k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive = FALSE, use.names = FALSE)))
}
