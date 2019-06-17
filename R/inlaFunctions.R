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
#' @param MRAcovParasGammaAlphaBeta list with two components, named 'space' and 'time'. Each component is itself a list with three components, 'rho', 'smoothness', and 'scale'. Finally, each of these thress components is a vector with two elements, corresponding to the alpha and beta parameters of the Gamma distribution serving as a prior for the MRA space and time hyperparameters.
#'
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

MRA_INLA <- function(spacetimeData, errorSDstart, fixedEffSDstart, MRAhyperparasStart, FEmuVec, predictionData = NULL, fixedEffGammaAlphaBeta, errorGammaAlphaBeta,  MRAcovParasGammaAlphaBeta, control) {
  defaultControl <- list(M = 1, randomSeed = 24,  cutForTimeSplit = 400, stepSize = 1, lowerThreshold = 3, maternCovariance = TRUE, nuggetSD = 0.00001, varyFixedEffSD = FALSE, varyMaternSmoothness = FALSE, varyErrorSD = TRUE, splitTime = FALSE, numKnotsRes0 = 20L, J = 2L, numValuesForGrid = 4)
  if (length(position <- grep(colnames(spacetimeData@sp@coords), pattern = "lon")) >= 1) {
    colnames(spacetimeData@sp@coords)[[position[[1]]]] <- "x"
  }
  if (length(position <- grep(colnames(spacetimeData@sp@coords), pattern = "lat")) >= 1) {
    colnames(spacetimeData@sp@coords)[[position[[1]]]] <- "y"
  }
  coordRanges <- mapply(dimName = c("lonRange", "latRange", "timeRange"), code = c("x", "y", "time"), function(dimName, code) {
    predCoordinates <- c()
    if (code != "time") {
      bufferSize <- 0.01
      coordinates <- spacetimeData@sp@coords[, code]
      if (!is.null(predictionData)) {
        predCoordinates <- predictionData@sp@coords[, code]
      }
    } else {
      bufferSize <- 10
      coordinates <- time(spacetimeData@time)
      if (!is.null(predictionData)) {
        predCoordinates <- time(predictionData@time)
      }
    }
    combinedRange <- range(c(coordinates, predCoordinates))
    combinedRange + c(-bufferSize, bufferSize)
  }, SIMPLIFY = FALSE)
  defaultControl <- c(defaultControl, coordRanges)

  for (i in names(control)) {
    defaultControl[[i]] <- control[[i]]
  }

  control <- defaultControl
  dataCoordinates <- spacetimeData@sp@coords
  timeRangeReshaped <- as.integer(control$timeRange)/(3600*24)
  timeBaseline <- min(timeRangeReshaped)
  timeValues <- as.integer(time(spacetimeData@time))/(3600*24) - timeBaseline # The division is to obtain values in days.
  timeRangeReshaped <- timeRangeReshaped - timeBaseline

  covariateMatrix <- as.matrix(spacetimeData@data[, -1, drop = FALSE])
  gridPointer <- setupGridCpp(spacetimeData@data[, 1], dataCoordinates, timeValues, covariateMatrix, control$M, control$lonRange, control$latRange, timeRangeReshaped, control$randomSeed, control$cutForTimeSplit, control$splitTime, control$numKnotsRes0, control$J)
  # First we compute values relating to the hyperprior marginal distribution...
  xStartValues <- c(spRho = MRAhyperparasStart$space[["rho"]], timeRho = MRAhyperparasStart$time[["rho"]], scale = MRAhyperparasStart[["scale"]])
  if (control$varyFixedEffSD) {
    xStartValues[["fixedEffSD"]] <- fixedEffSDstart
  }
  if (control$varyErrorSD) {
    xStartValues[["errorSD"]] <- errorSDstart
  }
  if (control$varyMaternSmoothness) {
    xStartValues[["spSmoothness"]] <- MRAhyperparasStart$space[["smoothness"]]
    xStartValues[["timeSmoothness"]] <- MRAhyperparasStart$time[["smoothness"]]
  }

  computedValues <- obtainGridValues(treePointer = gridPointer$gridPointer, xStartValues = xStartValues, control = control, fixedEffSDstart = fixedEffSDstart, errorSDstart = errorSDstart, MRAhyperparasStart = MRAhyperparasStart, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, FEmuVec = FEmuVec)

  discreteLogJointValues <- sapply(computedValues, '[[', "logJointValue")
  print(discreteLogJointValues)
  cat("\n \n")
  print(max(discreteLogJointValues))
  cat("\n \n")
  maxLogJointValues <- max(discreteLogJointValues)
  logStandardisingConstant <- maxLogJointValues + log(sum(exp(discreteLogJointValues - maxLogJointValues)))

  # The following normalises the joint distribution.
  cat("Standardising empirical distribution...\n") ;
  hyperparaList <- lapply(computedValues, function(x) {
    x$logJointValue <- x$logJointValue - logStandardisingConstant
    x
  })

  # Now, we obtain the marginal distribution of all mean parameters.
  cat("Computing moments for marginal posterior distributions...\n")
  hyperMarginalMoments <- ComputeHyperMarginalMoments(hyperparaList)
  meanMarginalMoments <- ComputeMeanMarginalMoments(hyperparaList)
  outputList <- list(hyperMarginalMoments = hyperMarginalMoments, meanMarginalMoments = meanMarginalMoments)
  cat("Computing prediction moments... \n")
  if (!is.null(predictionData)) {
    outputList$predictionMoments <- ComputeKrigingMoments(hyperparaList)
  }
  cat("Returning results... \n")
  outputList
}

obtainGridValues <- function(treePointer, xStartValues, control, fixedEffSDstart, errorSDstart, MRAhyperparasStart, MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, FEmuVec) {

  currentCenter <- xStartValues
  radius <- 0.9 * xStartValues
  currentMax <- -Inf
  containerList <- list()
  for (i in 1:20) {
    lowerLimit <- sapply(currentCenter - radius, function(x) ifelse(x <= 0.005, 0.005, x))
    upperLimit <- currentCenter + radius
    paraLimits <- cbind(lower = lowerLimit, upper = upperLimit)
    paraRanges <- lapply(1:nrow(paraLimits), function(x) seq(from = paraLimits[x,"lower"], to = paraLimits[x,"upper"], length.out = 2))
    names(paraRanges) <- names(currentCenter)
    paraGrid <- expand.grid(paraRanges)
    valuesOnGrid <- lapply(1:nrow(paraGrid), funForGridEst, paraGrid = paraGrid, treePointer = treePointer, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, fixedEffSDstart = fixedEffSDstart, errorSDstart = errorSDstart, MRAhyperparasStart = MRAhyperparasStart, FEmuVec = FEmuVec, control = control)
    keepIndices <- sapply(valuesOnGrid, function(x) class(x$logJointValue) == "numeric")
    valuesOnGrid <- valuesOnGrid[keepIndices]
    containerList <- c(containerList, valuesOnGrid)
    postProbValues <- sapply(valuesOnGrid, function(x) x$logJointValue)
    proposedMax <- max(postProbValues)

    if (proposedMax > currentMax) {
      currentCenter <- valuesOnGrid[[which.max(postProbValues)]]$x
      currentMax <- proposedMax
    } else {
      radius <- 0.75 * radius
      cat("Radius vector is now: ", radius, "\n")
    }
  }
  containerList
}

funForGridEst <- function(index, paraGrid, treePointer, MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, fixedEffSDstart, errorSDstart, MRAhyperparasStart, FEmuVec, control) {
  x <- unlist(paraGrid[index, ])
  fixedEffArg <- fixedEffSDstart
  if (control$varyFixedEffSD) {
    fixedEffArg <- x[["fixedEffSD"]]
  }
  errorArg <- errorSDstart
  if (control$varyErrorSD) {
    errorArg <- x[["errorSD"]]
  }
  spSmoothnessArg <- MRAhyperparasStart$space[["smoothness"]]
  timeSmoothnessArg <- MRAhyperparasStart$time[["smoothness"]]
  if (control$varyMaternSmoothness) {
    spSmoothnessArg <- x[["spSmoothness"]]
    timeSmoothnessArg <- x[["timeSmoothness"]]
  }
  MRAlist <- list(space = list(rho = x[["spRho"]], smoothness = spSmoothnessArg), time = list(rho = x[["timeRho"]], smoothness = timeSmoothnessArg), scale = x[["scale"]])
  logJointValue <- tryCatch(expr = LogJointHyperMarginal(treePointer = treePointer, MRAhyperparas = MRAlist, fixedEffSD = fixedEffArg, errorSD = errorArg, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = control$maternCovariance, spaceNuggetSD = control$nuggetSD, timeNuggetSD = control$nuggetSD, recordFullConditional = FALSE), error = function(e) e)
  errorSD <- ifelse(any(grepl(pattern = "error", x = names(x))), x[[grep(pattern = "error", x = names(x), value = TRUE)]], errorSDstart)
  fixedEffSD <- ifelse(any(grepl(pattern = "fixedEff", x = names(x))), x[[grep(pattern = "fixedEff", x = names(x), value = TRUE)]], fixedEffSDstart)
  aList <- list(x = x, errorSD = errorSD, fixedEffSD = fixedEffSD, MRAhyperparas = MRAlist, logJointValue = logJointValue)
  aList$FullCondMean <- GetFullCondMean(treePointer)
  aList$FullCondSDs <- GetFullCondSDs(treePointer)
  aList
}

ComputeHyperMarginalMoments <- function(hyperparaList) {
  domainCheck <- sapply(hyperparaList, function(x) x$logJointValue > -Inf)
  hyperparaList <- hyperparaList[domainCheck]
  psiAndMargDistMatrix <- t(sapply(hyperparaList, function(x) c(unlist(x$MRAhyperparas), fixedEffSD = x$fixedEffSD, errorSD = x$errorSD, jointValue = exp(x$logJointValue))))
  print(psiAndMargDistMatrix)
  rownames(psiAndMargDistMatrix) <- NULL
  funToGetParaMoments <- function(hyperparaIndex) {
    # variableValues <- sort(unique(psiAndMargDistMatrix[, hyperparaIndex]))
    # massValues <- do.call("c", by(psiAndMargDistMatrix, INDICES = variableValues, FUN = function(block) sum(block[, ncol(block)]), simplify = FALSE))
    meanValue <- sum(psiAndMargDistMatrix[, hyperparaIndex] * psiAndMargDistMatrix[, ncol(psiAndMargDistMatrix)])
    sdValue <- sqrt(sum(psiAndMargDistMatrix[, hyperparaIndex]^2 * psiAndMargDistMatrix[, ncol(psiAndMargDistMatrix)]) - meanValue^2)
    c(mean = meanValue, StdDev = sdValue)
  }
  paraMoments <- t(sapply(1:(ncol(psiAndMargDistMatrix) - 1), FUN = funToGetParaMoments))

  colnames(paraMoments) <- c("Mean", "StdDev")
  rownames(paraMoments) <- head(colnames(psiAndMargDistMatrix), n = -1)
  as.data.frame(paraMoments)
}

ComputeMeanMarginalMoments <- function(hyperparaList) {
  domainCheck <- sapply(hyperparaList, function(x) x$logJointValue > -Inf)
  hyperparaList <- hyperparaList[domainCheck]
  numMeanParas <- length(hyperparaList[[1]]$FullCondMean)
  weights <- exp(sapply(hyperparaList, "[[", "logJointValue"))
  marginalMeans <- sapply(1:numMeanParas, function(paraIndex) {
    meanVector <- sapply(hyperparaList, function(x) x$FullCondMean[[paraIndex]])
    sum(meanVector * weights)
  })
  marginalSecondMoments <- sapply(1:numMeanParas, function(paraIndex) {
    meanVector <- sapply(hyperparaList, function(x) x$FullCondMean[[paraIndex]])
    sdVector <- sapply(hyperparaList, function(x) x$FullCondSDs[[paraIndex]])
    secondMomentVec <- sdVector^2 + meanVector^2
    sum(secondMomentVec * weights)
  })
  marginalSDs <- marginalSecondMoments - marginalMeans^2
  data.frame(Mean = marginalMeans, StdDev = marginalSDs)
}

ComputeKrigingMoments <- function(hyperparaList) {
  domainCheck <- sapply(hyperparaList, function(x) x$logJointValue > -Inf)
  hyperparaList <- hyperparaList[domainCheck]
  termsForMean <- lapply(hyperparaList, function(x) {
    drop(x$CondPredStats$Hmean * exp(x$logJointValue))
  })
  krigingMeans <- Reduce("+", termsForMean)
  termsForVarE1 <- lapply(hyperparaList, FUN = function(x) {
    drop(x$CondPredStats$Hmean^2 * exp(x$logJointValue))
  })
  varE <- Reduce('+', termsForVarE1) - krigingMeans^2

  termsForEvar <- lapply(hyperparaList, function(x) {
    drop(x$CondPredStats$Evar * exp(x$logJointValue))
  })
  Evar <- Reduce("+", termsForEvar)
  list(predictMeans = krigingMeans, predictSDs = sqrt(varE + Evar))
}

covFunctionBiMatern <- function(rangeParaSpace = 10, rangeParaTime = 10) {
  function(spacetime1, spacetime2) {
    euclidDist <- sp::spDists(spacetime1@sp, spacetime2@sp)
    timeDist <- outer(zoo::index(spacetime2@time), zoo::index(spacetime1@time), function(x, y) as.numeric(abs(difftime(x, y, units = "days"))))
    fields::Exponential(euclidDist, range = 10)*t(fields::Exponential(timeDist, range = 10))
  }
}

# When lowerThreshold is 3, the exploration stops if the value being tested is at least 20 times smaller than the value at the peak of the hyperparameter marginal a posteriori distribution.
getIntegrationPointsAndValues <- function(optimObject, gridPointer, MRAcovParasGammaAlphaBeta, FEmuVec, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, stepSize = 1, lowerThreshold = 3, matern = FALSE, predictionData = NULL, timeBaseline, otherGridPointers = NULL, clusterAddress = NULL, spaceNuggetSD, timeNuggetSD, MRAhyperparasStart, errorSDstart, fixedEffSDstart) {

  decomposition <- eigen(solve(-optimObject$hessian), symmetric = TRUE)

  if (any(decomposition$values < 0)) {
    stop("Error: The negative Hessian matrix has negative eigenvalues. The values found in the optimisation is not a maximum. \n \n")
  }

  getPsi <- function(z) {
    sqrtEigenValueMatrix <- base::diag(sqrt(decomposition$values))
    drop(optimObject$par + decomposition$vectors %*% sqrtEigenValueMatrix %*% z)
  }

  getContainerElement <- function(z) {
    Psis <- getPsi(z)
    if (any(Psis <= 0)) {
      return(list(logJointValue = -Inf))
    }
    names(Psis) <- names(optimObject$par)
    MRAhyperList <- MRAhyperparasStart
    # spaceVars <- grep(pattern = "sp", x = names(z), value = TRUE)
    # timeVars <- grep(pattern = "time", x = names(z), value = TRUE)
    for (j in names(Psis)) {
      if (j == "spRho") {
        MRAhyperList$space[["rho"]] <- Psis[[j]]
      }
      if (j == "spSmoothness") {
        MRAhyperList$space[["smoothness"]] <- Psis[[j]]
      }
      if (j == "timeRho") {
        MRAhyperList$time[["rho"]] <- Psis[[j]]
      }
      if (j == "timeSmoothness") {
        MRAhyperList$time[["smoothness"]] <- Psis[[j]]
      }
      if (j == "scale") {
        MRAhyperList$scale <- Psis[[j]]
      }
    }
    errorSD <- ifelse(any(grepl(pattern = "error", x = names(Psis))), Psis[[grep(pattern = "error", x = names(Psis), value = TRUE)]], errorSDstart)
    fixedEffSD <- ifelse(any(grepl(pattern = "fixedEff", x = names(Psis))), Psis[[grep(pattern = "fixedEff", x = names(Psis), value = TRUE)]], fixedEffSDstart)

    aList <- list(z = z, errorSD = errorSD, fixedEffSD = fixedEffSD, MRAhyperparas = MRAhyperList)

    aList$logJointValue <- LogJointHyperMarginal(treePointer = gridPointer, MRAhyperparas = MRAhyperList, fixedEffSD = fixedEffSD, errorSD = errorSD, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta =  fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = matern, spaceNuggetSD = spaceNuggetSD, timeNuggetSD = timeNuggetSD, recordFullConditional = TRUE)
    if (!is.null(predictionData)) {
      timeValues <- as.integer(time(predictionData@time))/(3600*24) - timeBaseline # The division is to obtain values in days.
      aList$CondPredStats <- ComputeCondPredStats(gridPointer, predictionData@sp@coords, timeValues, as.matrix(predictionData@data))
    }
    # Running LogJointHyperMarginal stores in the tree pointed by gridPointer the full conditional mean and SDs when recordFullConditional = TRUE. We can get them with the simple functions I call now.
    aList$FullCondMean <- GetFullCondMean(gridPointer)
    aList$FullCondSDs <- GetFullCondSDs(gridPointer)
    aList
  }

  numDims <- length(optimObject$par)
  centerList <- list() # We'll be growing this list, but considering that we'll be having only a few hundred elements in final, this should be inconsequential.

  centerList[[1]] <- getContainerElement(rep(0,numDims))
  peakLogJointValue <- centerList[[1]]$logJointValue
  counter <- 1
  # This explores the distribution strictly along the axes defining centers for exploration at points not on the axes.
  exploreCombinations <- expand.grid(dimNumber = 1:numDims, direction = c(-1, 1))
  exploreList <- lapply(as.list(as.data.frame(t(exploreCombinations))), function(x) {
    names(x) <- colnames(exploreCombinations)
    x
  })
  # for (dimNumber in 1:numDims) {
  #   for (direction in c(-1, 1)) {
  centerListList <- lapply(exploreList, function(x) {
    centers <- list()
    counter <- 0
    currentZ <- rep(0, numDims)
    currentZ[x[["dimNumber"]]] <- stepSize*x[["direction"]]
    repeat {
      elementToAdd <- getContainerElement(currentZ)
      # The following check compares the value at the peak of the distribution with that at the current location. Since values are on the log-scale, the criterion is multiplicative: if the log of the ratio of the two density values is greater than the threshold, it means there was a steep enough drop, and that the density at the point considered is low enough that it won't matter much in the computation of the normalizing constant.
      if ((peakLogJointValue - elementToAdd$logJointValue) > lowerThreshold) break

      counter <- counter + 1
      centers[[counter]] <- elementToAdd
      currentZ[x[["dimNumber"]]] <- currentZ[x[["dimNumber"]]] + x[["direction"]]*stepSize
    }
    centers
  })
  #   }
  # }
  centerList <- do.call("c", centerListList)

  zMatrix <- t(sapply(centerList, '[[', "z"))
  containerList <- centerList
  counter <- length(containerList)
  # The following explores the distribution at points not found on the axes.
  for (centerIndex in 2:length(centerList)) {
    for (dimNumber in 1:numDims) {
      for (direction in c(-1, 1)) {
        currentZ <- centerList[[centerIndex]]$z
        repeat {
          currentZ[[dimNumber]] <- currentZ[[dimNumber]] + direction*stepSize
          # First, check if value is available
          rowNumber <- which(apply(zMatrix, MARGIN = 1, identical, currentZ))

          if (test <- length(rowNumber) == 0) { # This syntax actually works...
            containerElement <- getContainerElement(currentZ)
          } else {
            containerElement <- containerList[[rowNumber]]
          }

          if ((centerList[[1]]$logJointValue - containerElement$logJointValue) > lowerThreshold) break

          if (test) {
            counter <- counter + 1
            containerList[[counter]] <- containerElement
            zMatrix <- rbind(zMatrix, currentZ)
          }
        }
      }
    }
  }
  containerList
}

# SimulateSpacetimeData <- function(numObsPerTimeSlice = 2025, covFunction, lonRange, latRange, timeValuesInPOSIXct, covariateDistributionList, errorSD, distFun = dist, FEvalues, tiltSpaceSD = 0, tiltTime = FALSE) {
#   numSlotsPerRow <- ceiling(sqrt(numObsPerTimeSlice))
#   slotCoordinates <- sapply(list(longitude = lonRange, latitude = latRange), FUN = function(x) {
#     width <- abs(diff(x)/numSlotsPerRow)
#     seq(from = min(x) + width/2, to = max(x) - width/2, length.out = numSlotsPerRow)
#   }, USE.NAMES = TRUE)
#
#   allSpaceCoordinates <- as.data.frame(expand.grid(as.data.frame(slotCoordinates)))
#
#   timeLayerFct <- function(timeValue) {
#     selectedRows <- sample.int(n = nrow(allSpaceCoordinates), size = numObsPerTimeSlice, replace = FALSE)
#     layerCoordinates <- allSpaceCoordinates[selectedRows, ]
#     layerCoordinates$time <- rep(timeValue, nrow(layerCoordinates))
#     if (tiltTime) {
#       numObs <- length(selectedRows)
#       numObsDiv2 <- floor(numObs/2)
#       layerCoordinates$time <- layerCoordinates$time +  seq(from = -numObsDiv2, to = numObsDiv2 - ((numObs %% 2) == 0))
#     }
#     covariateFrame <- as.data.frame(sapply(covariateDistributionList, function(x) sample(x[[1]], size = nrow(layerCoordinates), prob = x[[2]], replace = TRUE)))
#     colnames(covariateFrame) <- paste("Covariate", seq_along(covariateFrame), sep = "")
#     cbind(layerCoordinates, covariateFrame)
#   }
#
#   coordinates <- lapply(timeValuesInPOSIXct, FUN = timeLayerFct)
#   coordinates <- do.call("rbind", coordinates)
#
#   spatialDistMatrix <- dist(coordinates[, c("longitude", "latitude")])
#   timeDistMatrix <- dist(coordinates[, "time"])/(3600*24)
#   covarianceMat <- covFunction(spatialDistMatrix, timeDistMatrix)
#   meanVector <- cbind(1, as.matrix(coordinates[ , grep(pattern = "Covariate", colnames(coordinates))])) %*% FEvalues
#   fieldValues <- mvtnorm::rmvnorm(n = 1, mean = as.numeric(meanVector), sigma = covarianceMat) + rnorm(n = length(meanVector), mean = 0, sd = errorSD)
#   dataForObject <- cbind(data.frame(y = fieldValues[1, ]), coordinates[ , grep(pattern = "Covariate", colnames(coordinates))])
#   spacetimeObj <- spacetime::STIDF(sp = sp::SpatialPoints(coordinates[, c("longitude", "latitude")]), data = dataForObject, time = coordinates$time)
#   if (tiltSpaceSD > 0) {
#     spacetimeObj@sp@coords <- spacetimeObj@sp@coords + rnorm(length(spacetimeObj@sp@coords), 0, tiltSpaceSD)
#   }
#   spacetimeObj
# }

SimulateSpacetimeData <- function(numObsPerTimeSlice = 225, covFunction, lonRange, latRange, timeValuesInPOSIXct, covariateGenerationFctList, errorSD, distFun = dist, FEvalues) {
  numSlotsPerRow <- ceiling(sqrt(numObsPerTimeSlice))
  slotCoordinates <- sapply(list(longitude = lonRange, latitude = latRange), FUN = function(x) {
    width <- abs(diff(x)/numSlotsPerRow)
    seq(from = min(x) + width/2, to = max(x) - width/2, length.out = numSlotsPerRow)
  }, USE.NAMES = TRUE)

  allSpaceCoordinates <- as.data.frame(expand.grid(as.data.frame(slotCoordinates)))
  numToRemove <- nrow(allSpaceCoordinates) - numObsPerTimeSlice
  obsToRemove <- (nrow(allSpaceCoordinates) - numToRemove + 1):nrow(allSpaceCoordinates)
  allSpaceCoordinates <- allSpaceCoordinates[-obsToRemove, ]

  coordinates <- allSpaceCoordinates[rep(1:nrow(allSpaceCoordinates), length(timeValuesInPOSIXct)), ]
  coordinates$time <- rep(timeValuesInPOSIXct, each = numObsPerTimeSlice)

  covariateMatrix <- cbind(1, as.matrix(sapply(covariateGenerationFctList, function(x) x(coordinates[, c("longitude", "latitude")], coordinates[, "time"]))))

  spatialDistMatrix <- dist(coordinates[, c("longitude", "latitude")])
  timeDistMatrix <- dist(coordinates[, "time"])/(3600*24)
  covarianceMat <- covFunction(spatialDistMatrix, timeDistMatrix)
  meanVector <- drop(covariateMatrix %*% FEvalues)
  fieldValues <- drop(mvtnorm::rmvnorm(n = 1, mean = meanVector, sigma = covarianceMat)) + rnorm(n = length(meanVector), mean = 0, sd = errorSD)
  dataForObject <- cbind(y = fieldValues, as.data.frame(covariateMatrix[, -1]))
  colnames(dataForObject) <- c("y", paste("Covariate", 1:(length(FEvalues) - 1), sep = ""))
  spacetimeObj <- spacetime::STIDF(sp = sp::SpatialPoints(coordinates[, c("longitude", "latitude")]), data = dataForObject, time = coordinates$time)
  spacetimeObj
}

# In the Wikipedia notation, smoothness corresponds to nu, and
# scale corresponds to sigma.
maternCov <- function(d, rho, smoothness, scale) {
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d[d == 0] <- 1e-10

  dScaled <- sqrt(2 * smoothness) * d / rho
  con <- scale^2 * 2^(1 - smoothness) / gamma(smoothness)

  con * dScaled^smoothness * besselK(dScaled, smoothness)
}

LogJointHyperMarginal <- function(treePointer, MRAhyperparas, fixedEffSD, errorSD, MRAcovParasGammaAlphaBeta, FEmuVec, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, matern, spaceNuggetSD, timeNuggetSD, recordFullConditional) {
  # Hmat is the covariate matrix with a column of 1s at the front for the intercept, with a n x n identity matrix horizontally appended (horizontal/row merge).
  # MRAprecision has to be a sparse matrix.
  require(Matrix, quietly = TRUE)
  choleskiSolve <- function(hessianMat, scaledResponseHmat) {
    hessianMat <- as(hessianMat, "symmetricMatrix")
    value <- as.vector(Matrix::solve(hessianMat, as.vector(scaledResponseHmat), sparse = TRUE))
    value
  }

  funForTrustOptim <- function(start, sigmaSqEpsilon, MRAprecision, responseVec, Hmat) {
    cat("Entered funForTrustOptim...")

    funForOptim <- function(x) {
      recenteredResponse <- responseVec - Hmat %*% x
      firstTerm <- t(x) %*% MRAprecision %*% x
      secondTerm <- 1/sigmaSqEpsilon * Matrix::t(recenteredResponse) %*% recenteredResponse
      value <- -0.5 * (firstTerm + secondTerm)
      value[1, 1]
    }

    hessianMat <- -MRAprecision - 1/sigmaSqEpsilon * Matrix::t(Hmat) %*% Hmat
    intercept <- 1/sigmaSqEpsilon * t(responseVec) %*% Hmat

    gradForOptim <- function(x) {
      value <- t(x) %*% hessianMat + intercept
      value[1, ]
    }

    if (!("dgCMatrix" %in% class(hessianMat))) {
      stop("R error: Hessian is supposed to be sparse. \n")
    }

    hessianForOptim <- function(x) {
      hessianMat
    }

    opt <- trustOptim::trust.optim(x = drop(start), fn = funForOptim,
                       gr = gradForOptim,
                       hs = hessianForOptim,
                       method = "Sparse",
                       control = list(
                         start.trust.radius = 3,
                         stop.trust.radius = 1e-4,
                         prec = 1e-4,
                         report.precision = 1L,
                         maxit = 15L, # Since we are solving a linear problem, convergence is in one iteration.
                         preconditioner = 0,
                         precond.refresh.freq = 1,
                         function.scale.factor = -1 # We are maximising, hence the -1.
                       ))
    cat("Partial solution: ", opt$solution[1:8])
    opt$solution
  }
  LogJointHyperMarginalToWrap(treePointer = treePointer, MRAhyperparas = MRAhyperparas, fixedEffSD = fixedEffSD, errorSD = errorSD, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta =  fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = matern, spaceNuggetSD = spaceNuggetSD, timeNuggetSD = timeNuggetSD, recordFullConditional = TRUE, optimFun = funForTrustOptim, gradCholeskiFun = choleskiSolve)
}
