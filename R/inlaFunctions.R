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

MRA_INLA <- function(spacetimeData, errorSDstart, fixedEffSDstart, MRAhyperparasStart, M, lonRange, latRange, timeRange, randomSeed, cutForTimeSplit = 400, MRAcovParasIGalphaBeta, FEmuVec, fixedEffIGalphaBeta, errorIGalphaBeta, stepSize = 1, lowerThreshold = 10) {
  dataCoordinates <- spacetimeData@sp@coords
  timeValues <- as.integer(time(spacetimeData@time))/(3600*24) # The division is to obtain values in days.
  covariateMatrix <- as.matrix(spacetimeData@data[, -1, drop = FALSE])
  gridPointer <- setupGridCpp(spacetimeData@data[, 1], dataCoordinates, timeValues, covariateMatrix, M, lonRange, latRange, as.integer(timeRange)/(3600*24) + 1, randomSeed, cutForTimeSplit)
  # First we compute values relating to the hyperprior marginal distribution...

  xStartValues <- c(MRAhyperparasStart, fixedEffSDstart, errorSDstart)
  numMRAhyperparas <- length(MRAhyperparasStart)
  # funForOptimJointHyperMarginal(treePointer = gridPointer$gridPointer, exp(0.5*xStartValues[1:numMRAhyperparas]), exp(0.5*xStartValues[numMRAhyperparas + 1]), exp(0.5*xStartValues[numMRAhyperparas + 2]), MRAcovParasIGalphaBeta = MRAcovParasIGalphaBeta, fixedEffIGalphaBeta = fixedEffIGalphaBeta, errorIGalphaBeta = errorIGalphaBeta)
  funForOptim <- function(x, treePointer, MRAcovParasIGalphaBeta, fixedEffIGalphaBeta, errorIGalphaBeta, FEmuVec) {
    -LogJointHyperMarginal(treePointer, x[1:numMRAhyperparas], x[numMRAhyperparas + 1], x[numMRAhyperparas + 2], MRAcovParasIGalphaBeta, FEmuVec, fixedEffIGalphaBeta, errorIGalphaBeta)
  }
  optimResult <- optim(par = xStartValues, method = "L-BFGS-B", control = list(factr = 1e6), hessian = TRUE, fn = funForOptim, lower = rep(0.1, 4), gr = NULL, treePointer = gridPointer$gridPointer, MRAcovParasIGalphaBeta = MRAcovParasIGalphaBeta, FEmuVec = FEmuVec, fixedEffIGalphaBeta = fixedEffIGalphaBeta, errorIGalphaBeta = errorIGalphaBeta)
  cat("Best vector: ", optimResult$par, "\n")
  cat("Optimal value: ", optimResult$value, "\n")
  if (optimResult$convergence > 0) {
    warning("Nelder-Mead failed to conclusively find a maximum for the marginal hyperpara. posterior. \n Algorithm will keep going using the best value found. \n \n")
  }
  # save(optimResult, file = "~/Documents/optimForTests.R")
  # load("~/Documents/optimForTests.R")
  # cat("LOADING VALUES FOR TESTING PURPOSES. REMOVE THIS ONCE CODE IS COMPLETE.\n \n")
  # hyperparaList <- getIntegrationPointsAndValues(optimObject = optimResult, gridPointer = gridPointer$gridPointer, MRAcovParasIGalphaBeta = MRAcovParasIGalphaBeta, fixedEffIGalphaBeta = fixedEffIGalphaBeta, errorIGalphaBeta = errorIGalphaBeta, stepSize = stepSize, lowerThreshold = lowerThreshold)
  # standardisingConstant <- sum(sapply(hyperparaList, `$`, "logJointValue"))
  # # The following normalises the joint distribution.
  # hyperparaList <- lapply(hyperparaList, function(x) {
  #   x$logJointValue <- x$logJointValue - standardisingConstant
  #   x
  # })
  # # Now, we obtain the marginal distribution of all mean parameters.
  # ComputeMarginals(hyperparaList)
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

covFunctionBiMatern <- function(rangeParaSpace = 10, rangeParaTime = 10) {
  function(spacetime1, spacetime2) {
    euclidDist <- spDists(spacetime1@sp, spacetime2@sp)
    timeDist <- outer(zoo::index(spacetime2@time), zoo::index(spacetime1@time), function(x, y) as.numeric(abs(difftime(x, y, units = "days"))))
    fields::Exponential(euclidDist, range = 10)*t(fields::Exponential(timeDist, range = 10))
  }
}

getIntegrationPointsAndValues <- function(optimObject, gridPointer, MRAcovParasIGalphaBeta, fixedEffIGalphaBeta, errorIGalphaBeta, stepSize = 1, lowerThreshold = 10) {

  decomposition <- eigen(solve(-optimObject$hessian), symmetric = TRUE)

  getPsi <- function(z) {
    sqrtEigenValueMatrix <- diag(sqrt(decomposition$values))
    optimObject$value + decomposition$vectors %*% sqrtEigenValueMatrix %*% z
  }

  getContainerElement <- function(z) {
    Psis <- getPsi(z) ;
    aList <- vector(mode = 'list', length = 1)
    aList$z <- z
    aList$MRAhyperparas <- head(Psis, n = -2)
    aList$fixedEffSD <- Psis[length(Psis) - 1]
    aList$errorSD <- tail(Psis, n = 1)
    aList$logJointValue <- with(aList, expr = LogJointHyperMarginal(gridPointer, MRAhyperparas, fixedEffSD = fixedEffSD, errorSD = errorSD, MRAcovParasIGalphaBeta, fixedEffIGalphaBeta, errorIGalphaBeta))
    aList
  }

  numDims <- length(optimObject$par)
  centerList <- vector(mode = 'list', length = 5000) # We should not need more than a few hundred points, so 5000 should be ok.

  centerList[[1]] <- getContainerElement(rep(0,numDims))
  counter <- 2
  previousValue <- centerList[[1]]
  for (dimNumber in 1:numDims) {
    for (direction in c(-1, 1)) {
      currentZ <- rep(0, numDims)
      currentZ[dimNumber] <- stepSize*direction
      repeat {
        centerList[[counter]] <- getContainerElement(currentZ)
        counter = counter + 1
        if (abs(centerList[[counter]]$jointValue - previousValue) < lowerThreshold) break

        previousValue <- centerList[[counter]]$jointValue
        currentZ[dimNumber] = currentZ[dimNumber] + direction*stepSize ;
      }
    }
  }
  zMatrix <- t(sapply(centerList, `$`, "z"))
  containerList <- centerList

  for (centerIndex in 2:length(centerList)) {
    for (dimNumber in 1:numDims) {
      for (direction in c(-1, 1)) {
        currentZ <- centerList[[centerIndex]]$z
        previousValue <- centerList[[centerIndex]]$logJointValue
        repeat {
          currentZ[[dimNumber]] <- currentZ[[dimNumber]] + direction*stepSize
          # First, check if value is available
          rowNumber <- which(apply(zMatrix, MARGIN = 1, all.equal, currentZ))
          if (length(rowNumber) == 0) {
            containerElement <- getContainerElement(currentZ)
            containerList[[counter]] <- containerElement
            zMatrix <- rbind(zMatrix, currentZ)
            counter <- counter + 1
          } else {
            containerElement <- containerList[[rowNumber]]$logJointValue
          }
          if (abs(containerElement - previousValue) < lowerThreshold) break
          previousValue <- containerElement
        }
      }
    }
  }
  containerList[1:counter]
}

ComputeMarginals <- function(hyperparameterDistList) {
  dfsfsd
}

SimulateSpacetimeData <- function(numObsPerTimeSlice = 2025, covFunction, lonRange, latRange, timeValuesInPOSIXct, covariateDistributionList, errorSD, distFun = dist, FEvalues) {
  numSlotsPerRow <- ceiling(sqrt(numObsPerTimeSlice))
  slotCoordinates <- sapply(list(longitude = lonRange, latitude = latRange), FUN = function(x) {
    width <- abs(diff(x)/numSlotsPerRow)
    seq(from = min(x) + width/2, to = max(x) - width/2, length.out = numSlotsPerRow)
  }, USE.NAMES = TRUE)

  allSpaceCoordinates <- as.data.frame(expand.grid(as.data.frame(slotCoordinates)))

  timeLayerFct <- function(timeValue) {
    selectedRows <- sample.int(n = nrow(allSpaceCoordinates), size = numObsPerTimeSlice, replace = FALSE)
    layerCoordinates <- allSpaceCoordinates[selectedRows, ]
    layerCoordinates$time <- rep(timeValue, nrow(layerCoordinates))
    covariateFrame <- as.data.frame(sapply(covariateDistributionList, function(x) sample(x[[1]], size = nrow(layerCoordinates), prob = x[[2]], replace = TRUE)))
    colnames(covariateFrame) <- paste("Covariate", seq_along(covariateFrame), sep = "")
    cbind(layerCoordinates, covariateFrame)
  }

  coordinates <- lapply(timeValuesInPOSIXct, FUN = timeLayerFct)
  coordinates <- do.call("rbind", coordinates)

  spatialDistMatrix <- dist(coordinates[, c("longitude", "latitude")])
  timeDistMatrix <- dist(coordinates[, "time"])/(3600*24) # Distances in days
  covarianceMat <- covFunction(spatialDistMatrix, timeDistMatrix)
  meanVector <- cbind(1, as.matrix(coordinates[ , grep(pattern = "Covariate", colnames(coordinates))])) %*% FEvalues
  fieldValues <- mvtnorm::rmvnorm(n = 1, mean = as.numeric(meanVector), sigma = covarianceMat) + rnorm(n = length(meanVector), mean = 0, sd = errorSD)
  dataForObject <- cbind(data.frame(y = fieldValues[1, ]), coordinates[ , grep(pattern = "Covariate", colnames(coordinates))])
  spacetime::STIDF(sp = sp::SpatialPoints(coordinates[, c("longitude", "latitude")]), data = dataForObject, time = coordinates$time)
}
