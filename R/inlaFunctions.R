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

MRA_INLA <- function(spacetimeData, errorSDstart, fixedEffSDstart, MRAhyperparasStart, M, lonRange, latRange, timeRange, randomSeed, cutForTimeSplit = 400, MRAcovParasGammaAlphaBeta, FEmuVec, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, stepSize = 1, lowerThreshold = 3, maternCov = FALSE, spaceNuggetSD, timeNuggetSD, predictionData = NULL, clusterAddress = NULL, varyFixedEffSD = TRUE, varyMaternSmoothness = TRUE, varyErrorSD = TRUE, splitTime = TRUE) {
  dataCoordinates <- spacetimeData@sp@coords
  timeRangeReshaped <- as.integer(timeRange)/(3600*24)
  timeBaseline <- min(timeRangeReshaped)
  timeValues <- as.integer(time(spacetimeData@time))/(3600*24) - timeBaseline # The division is to obtain values in days.
  timeRangeReshaped <- timeRangeReshaped - timeBaseline

  covariateMatrix <- as.matrix(spacetimeData@data[, -1, drop = FALSE])
  gridPointer <- setupGridCpp(spacetimeData@data[, 1], dataCoordinates, timeValues, covariateMatrix, M, lonRange, latRange, timeRangeReshaped, randomSeed, cutForTimeSplit, splitTime)
  # First we compute values relating to the hyperprior marginal distribution...
  xStartValues <- c(spRho = MRAhyperparasStart$space[["rho"]], spScale = MRAhyperparasStart$space[["scale"]], timeRho = MRAhyperparasStart$time[["rho"]], timeScale = MRAhyperparasStart$time[["scale"]])
  if (varyFixedEffSD) {
    xStartValues[["fixedEffSD"]] <- fixedEffSDstart
  }
  if (varyErrorSD) {
    xStartValues[["errorSD"]] <- errorSDstart
  }
  if (varyMaternSmoothness) {
    xStartValues[["spSmoothness"]] <- MRAhyperparasStart$space[["smoothness"]]
    xStartValues[["timeSmoothness"]] <- MRAhyperparasStart$time[["smoothness"]]
  }

  funForOptim <- function(x, treePointer, MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, FEmuVec, elementNames) {
    names(x) <- elementNames
    fixedEffArg <- fixedEffSDstart
    if (varyFixedEffSD) {
      fixedEffArg <- x[["fixedEffSD"]]
    }
    errorArg <- errorSDstart
    if (varyErrorSD) {
      errorArg <- x[["errorSD"]]
    }
    spSmoothnessArg <- MRAhyperparasStart$space[["smoothness"]]
    timeSmoothnessArg <- MRAhyperparasStart$time[["smoothness"]]
    if (varyMaternSmoothness) {
      spSmoothnessArg <- x[["spSmoothness"]]
      timeSmoothnessArg <- x[["timeSmoothness"]]
    }
    MRAlist <- list(space = list(rho = x[["spRho"]], smoothness = spSmoothnessArg, scale = x[["spScale"]]), time = list(rho = x[["timeRho"]], smoothness = timeSmoothnessArg, scale = x[["timeScale"]]))
    -LogJointHyperMarginal(treePointer = treePointer, MRAhyperparas = MRAlist, fixedEffSD = fixedEffArg, errorSD = errorArg, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = as.logical(maternCov), spaceNuggetSD = spaceNuggetSD, timeNuggetSD = timeNuggetSD, recordFullConditional = FALSE)
  }

  cat("Optimising marginal hyperparameter posterior distribution... \n") ;
  optimResult <- nloptr::lbfgs(x0 = xStartValues, fn = funForOptim, lower = rep(0.01, length(xStartValues)), upper = 10*xStartValues, control = list(maxeval = 200, xtol_rel = 1e-4, ftol_rel = 1e-4), treePointer = gridPointer$gridPointer, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, elementNames = names(xStartValues), fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta)
  cat("Optimisation complete... \n")
  if (optimResult$convergence < 0) {
    stop("Optimisation algorithm did not converge! \n \n")
  }
  # save(optimResult, file = "~/Documents/optimForTests.R")
  # cat("LOADING VALUES FOR TESTING PURPOSES. REMOVE THIS ONCE CODE IS COMPLETE.\n \n")
  # load("~/Documents/optimForTests.R")
  cat("Estimating Hessian at mode... \n") ;
  hessianMat <- nlme::fdHess(pars = optimResult$par, fun = funForOptim, treePointer = gridPointer$gridPointer, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, elementNames = names(xStartValues), fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, .relStep = .Machine$double.eps^(1/3))$Hessian

  if (det(-hessianMat) <= 0) {
    warning("Non-positive definite negative Hessian matrix... \n \n")
  }
  optimResult$hessian <- -hessianMat # The "-" is necessary because we performed a minimisation rather than a maximisation.
  return(optimResult)
  cat("Computing distribution value at integration points... \n")
  hyperparaList <- getIntegrationPointsAndValues(optimObject = optimResult, gridPointer = gridPointer$gridPointer, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, stepSize = stepSize, lowerThreshold = lowerThreshold, matern = maternCov, predictionData = predictionData, timeBaseline = timeBaseline)
  # load("~/Documents/hyperparaList.Rdata")
  discreteLogJointValues <- sapply(hyperparaList, '[[', "logJointValue")
  maxLogJointValues <- max(discreteLogJointValues)
  logStandardisingConstant <- maxLogJointValues + log(sum(exp(discreteLogJointValues - maxLogJointValues)))

  # The following normalises the joint distribution.
  cat("Standardising empirical distribution...\n") ;
  hyperparaList <- lapply(hyperparaList, function(x) {
    x$logJointValue <- x$logJointValue - logStandardisingConstant
    x
  })

  # Now, we obtain the marginal distribution of all mean parameters.
  cat("Computing moments for marginal posterior distributions...\n")
  hyperMarginalMoments <- ComputeHyperMarginalMoments(hyperparaList)
  meanMarginalMoments <- ComputeMeanMarginalMoments(hyperparaList)
  outputList <- list(hyperMarginalMoments = hyperMarginalMoments, meanMarginalMoments = meanMarginalMoments)
  if (!is.null(predictionData)) {
    outputList$predictionMoments <- ComputeKrigingMoments(hyperparaList)
  }
  cat("Returning results... \n")
  outputList
}

ComputeHyperMarginalMoments <- function(hyperparaList) {
  psiAndMargDistMatrix <- t(sapply(hyperparaList, function(x) c(x$MRAhyperparas, x$fixedEffSD, x$errorSD, exp(x$logJointValue))))
  funToGetParaMoments <- function(hyperparaIndex) {
    variableValues <- sort(unique(psiAndMargDistMatrix[, hyperparaIndex]))
    massValues <- do.call("c", by(psiAndMargDistMatrix, INDICES = variableValues, FUN = function(block) sum(block[, ncol(block)]), simplify = FALSE))
    meanValue <- sum(variableValues * massValues)
    sdValue <- sqrt(sum(variableValues^2 * massValues) - meanValue^2)
    c(mean = meanValue, StdDev = sdValue)
  }
  paraMoments <- t(sapply(1:(ncol(psiAndMargDistMatrix) - 1), FUN = funToGetParaMoments))

  colnames(paraMoments) <- c("Mean", "StdDev")
  as.data.frame(paraMoments)
}

ComputeMeanMarginalMoments <- function(hyperparaList) {
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
getIntegrationPointsAndValues <- function(optimObject, gridPointer, MRAcovParasGammaAlphaBeta, FEmuVec, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, stepSize = 1, lowerThreshold = 3, matern = FALSE, predictionData = NULL, timeBaseline, otherGridPointers = NULL, clusterAddress = NULL) {

  decomposition <- eigen(solve(-optimObject$hessian), symmetric = TRUE)

  if (any(decomposition$values < 0)) {
    stop("Error: The Hessian matrix has negative eigenvalues. The values found in the optimisation is not a maximum. \n \n")
  }

  getPsi <- function(z) {
    sqrtEigenValueMatrix <- base::diag(sqrt(decomposition$values))
    drop(optimObject$par + decomposition$vectors %*% sqrtEigenValueMatrix %*% z)
  }

  getContainerElement <- function(z) {
    Psis <- getPsi(z)
    aList <- list(z = z, MRAhyperparas = head(Psis, n = -2), fixedEffSD = Psis[length(Psis) - 1], errorSD = tail(Psis, n = 1))
    aList$logJointValue <- with(aList, expr = LogJointHyperMarginal(treePointer = gridPointer, MRAhyperparas = MRAhyperparas, fixedEffSD = fixedEffSD, errorSD = errorSD, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta =  fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = matern, spaceNuggetSD = spaceNuggetSD, timeNuggetSD = timeNuggetSD, recordFullConditional = TRUE))
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
  centerList <- list() ; # We'll be growing this list, but considering that we'll be having only a few hundred elements in final, this should be inconsequential.

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
            cat("Point known... \n")
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
