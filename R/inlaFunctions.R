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

MRA_INLA <- function(spacetimeData, errorSDstart, fixedEffSDstart, MRAhyperparasStart, M, lonRange, latRange, timeRange, randomSeed, cutForTimeSplit = 400, hyperAlpha, hyperBeta) {
  dataCoordinates <- spacetimeData@sp@coords
  timeValues <- as.integer(time(spacetimeData@time))
  covariateMatrix <- as.matrix(spacetimeData@data[, -1, drop = FALSE])
  gridPointer <- setupGridCpp(spacetimeData@data[, 1], dataCoordinates, timeValues, covariateMatrix, M, lonRange, latRange, timeRange, randomSeed, cutForTimeSplit)
  # First we compute values relating to the hyperprior marginal distribution...
  # xStartValues <- log(c(MRAhyperparasStart, fixedEffSDstart, errorSDstart))
  # numMRAhyperparas <- length(MRAhyperparasStart)
  # funForOptim <- function(x, treePointer, hyperAlpha, hyperBeta) {
  #   result <- funForOptimJointHyperMarginal(treePointer, exp(x[1:numMRAhyperparas]), exp(x[numMRAhyperparas + 1]), exp(x[numMRAhyperparas + 2]), hyperAlpha, hyperBeta)
  #   result
  # }
  # valeur1 <- funForOptimJointHyperMarginal(gridPointer$gridPointer, MRAhyperparasStart,
  #                               fixedEffSDstart, errorSDstart, hyperAlpha, hyperBeta)
  # cat("Obtained valeur1! \n") ;
  valeur2 <- funForOptimJointHyperMarginal(gridPointer$gridPointer, MRAhyperparasStart,
                                           fixedEffSDstart, errorSDstart, hyperAlpha, hyperBeta)
  cat("Valeurs 2: ", valeur2, "\n")
  # cat("Start values: ", xStartValues, "\n")
  # optimResult <- optim(par = xStartValues, fn = funForOptim, gr = NULL, treePointer = gridPointer$gridPointer, hyperAlpha = hyperAlpha, hyperBeta = hyperBeta)
  # optimResult
  # hyperMode <- optimResult$par
  # hyperHessian <- optimResult$hessian
  # list(hyperDistMode = hyperMode, hessianAtMode = hyperHessian)
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
