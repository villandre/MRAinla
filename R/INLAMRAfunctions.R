#' INLA-MRA model for inference and prediction in spatiotemporal data
#'
#' The function fits the INLA-MRA model to a spatiotemporal dataset and outputs posterior predictive distributions. INLA-MRA assumes a multiplicative form for spatiotemporal covariance, with each component expressed with the Matérn formula.
#'
#' @param responseVec A numeric vector with response values.
#' @param covariateFrame A data.frame containing covariate values in the order of elements in responseVec
#' @param spatialCoordMat A matrix or data.frame with two columns with the *first corresponding to longitude*, and the *second to latitude*; can also be x and y if the sinusoidal projection is used (default for satellite imagery data)
#' @param timePOSIXctVec A vector of time values in POSIXct format
#' @param predCovariateFrame A data.frame containing covariate values for the prediction datasets
#' @param predSpatialCoordMat Like spatialCoordMat, but for the prediction data
#' @param predTimePOSIXctVec Like timePOSIXctVec, but for prediction data
#' @param sinusoidalProjection Logical value indicating whether the provided coordinates are in the sinusoidal projection
#' @param spatialRangeList List with two elements: a starting value for the *spatial range* hyperparameter, and a two element vector giving the mean and standard deviation of the normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param spatialSmoothnessList List with two elements: a starting value for the *spatial smoothness* hyperparameter, and a length-two vector giving the mean and standard deviation of the associated normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param timeRangeList List with two elements: a starting value for the *temporal range* hyperparameter, and a length-two vector giving the mean and standard deviation of the associated normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param timeSmoothnessList List with two elements: a starting value for the *temporal smoothness* hyperparameter, and a length-two vector giving the mean and standard deviation of the associated normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param scaleList List with two elements: a starting value for the *scale* hyperparameter, and a length-two vector giving the mean and standard deviation of the associated normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param errorSDlist List with two elements: a starting value for the *uncorrelated error standard deviation* hyperparameter, and a length-two vector giving the mean and standard deviation of the associated normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param fixedEffSDlist List with two elements: a starting value for the *uncorrelated fixed effects standard deviation* hyperparameter, and a length-two vector giving the mean and standard deviation of the associated normal hyperprior (second element must be omitted if hyperparameter is fixed)
#' @param FEmuVec Vector with the mean value of the priors for the fixed effects. Its length must match the number of columns in covariateFrame
#' @param control List with control parameters. See details.
#'
#' @details Some of the control parameters should be tuned to ensure better computational or predictive performance, or to make it possible to stop and resume model fitting:
#' \itemize{
#' \item{Mlon, Mlat, Mtime} {The number of longitude, latitude, and time splits used to create the nested resolutions. We have M = Mlon + Mlat + Mtime. They should be set as low as possible,  keeping in mind time and memory constraints.}
#' \item{numKnotsRes0} {The number of knots at resolution 0 (the resolution encompassing the entire spatiotemporal domain). Takes value 20 by default. We would not recommend setting it under 8, as knots are first placed on the vertices of a rectangular prism nested within each subregion. Increasing this number also increases the program's memory footprint and running time.}
#' \item{J} {Multiplier used to determine the number of knots at each resolution. The number of knots in each subregion at resolution i is ceiling(numKnotsRes0 * (J/2)^i). {Takes value 2 by default. We do not recommend setting J under 2, as some regions might end up with only two knots when M is large enough.}
#' \item{numValuesForIS} {The number of IS samples to be used for marginalisation. Takes value 100 by default.}
#' \item{numIterOptim} {The number of iterations in the L-BFGS algorithm used to identify the maximum of the join marginal posterior distribution of the hyperparameters. Takes value 25 by default. Could be set somewhat lower, 20 say, to reduce running time, but setting it too low might create imbalance in the importance sampling weights.}
#' \item{tipKnotsThinningRate} {The proportion of observation spatiotemporal locations that should be used as knots at the finest resolution, values should be in (0, 1]. Takes value 1 by default. A lower value for this parameter could reduce predictive performance, but greatly reduce the memory footprint. For very large datasets, lower values of this parameter are recommended, and can even be necessary.}
#' \item{credIntervalPercs} {The quantile boundaries of the reported credibility intervals. By default, 0.025 and 0.975, for a 95% credibility interval.}
#' \item{fileToSaveOptOutput} {String indicating where the results of the optimisation, used to identify the maximum of the marginal joint hyperparameter posterior distribution, should be saved. If this is set, it becomes possible to resume the fitting after interruption. The function will restart after the optimisation step.}
#' }
#' \item{folderToSaveISpoints} {String indicating the name of a folder where the results of each iteration of the IS algorithm should be saved. This allows the IS algorithm to be resumed after interruption. It also makes it possible to produce an output before all IS iterations have been processed, cf. control$IScompleted.}
#' There are other control parameters that should not required to be changed (but that we keep there in case they might be required in future versions of the package):
#' \itemize{
#' \item{distMethod} {String indicating the method used to obtain distances in kilometers from longitude/latitude coordinates. Takes value "haversine" by default, for the Haversine distance formula. No alternative is implemented for now.}
#' \item {randomSeed} {Seed for the random number generator used for knot placement. Takes value 24 by default.}
#' \item{numISpropDistUpdates} {The number of importance sampling (IS) weight updates in the adaptive IS algorithm. Takes value 0 by default. We implemented an adaptive IS scheme in case a reasonable level of balance in IS weights was not reached. It is enabled by setting this pararameter to 1 or more.}
#' \item{nuggetSD} {A small number added to the diagonal of the covariance matrices obtained by applying the Matern formula, to ensure that they are invertible. Takes value 1e-5 by default.}
#' \item{normalHyperprior} {Should hyperparameters be modelled on the logarithmic scale and normal hyperpriors be used? Takes value TRUE by default. Modelling hyperparameters on the original scale, with gamma priors, is possible, but not recommended.}
#' \item{IScompleted} {Logical value indicating whether all importance sampling weights have been obtained. Takes value FALSE by default. Set it to TRUE to produce intermediate results (based on fewer IS iterations than had been originally planned) when the run takes too long time to finish. Note that control$fileToSaveOptOutput and control$folderToSaveISpoints must be specified for this feature to work.}
#' \item{spaceJitterMax} {The maximum jittering to apply to longitude/ latitude coordinates. Takes value 0 by default, for no jittering. The method might run into numerical difficulties if longitude/latitude coordinates are replicated. The jittering ensures that it does not happen.}
#' }
#'
#'
#'
#' @return A list with three components:
#' \itemize{
#'  \item{hyperMarginalMoments} {A data.frame giving the mean, and standard deviation of the marginal hyperparameter posteriors, as well as the related 95\% credibility intervals.}
#'  \item{FEmarginalMoments} {A data.frame giving the mean, and standard deviation of the marginal fixed effects posteriors, as well as 95\% credibility intervals.}
#'  \item{predictionMoments} {A data.frame with two columns, predictMeans and predictSDs. The order of the predictions matches the one in predCovariateFrame.}
#' }
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

MRA_INLA <- function(responseVec, covariateFrame, spatialCoordMat, timePOSIXctVec, predCovariateFrame = NULL, predSpatialCoordMat = NULL, predTimePOSIXctVec = NULL, sinusoidalProjection = FALSE,  spatialRangeList = NULL, spatialSmoothnessList = list(start = log(1.5)), timeRangeList = NULL, timeSmoothnessList = list(start = log(0.5)), scaleList = list(start = 0, hyperpars = c(0, 2)), errorSDlist = list(start = 0), fixedEffSDlist = list(start = log(10)), FEmuVec = rep(0, ncol(covariateFrame)), control) {

  noPredictionFlag <- is.null(predCovariateFrame) | is.null(predSpatialCoordMat) | is.null(predTimePOSIXctVec)

  .checkInputConsistency(responseVec, covariateFrame, spatialCoordMat, timePOSIXctVec, predCovariateFrame, predSpatialCoordMat, predTimePOSIXctVec)

  # .checkInputConsistency has already ensured that column names match. This line ensures that covariates are presented in the exact same order in predictions as in observations.
  if (!noPredictionFlag) predCovariateFrame <- predCovariateFrame[colnames(covariateFrame)]
  ##################################
  lonLatProjString <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  crsString <- ifelse(sinusoidalProjection, yes = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs", no = lonLatProjString)
  spObject <- sp::SpatialPoints(coords = spatialCoordMat, proj4string = sp::CRS(crsString))
  spObjectPred <- NULL
  if (!noPredictionFlag) spObjectPred <- sp::SpatialPoints(coords = predSpatialCoordMat, proj4string = sp::CRS(crsString))
  if (sinusoidalProjection) {
    spObject <- sp::spTransform(x = spObject, CRSobj = lonLatProjString)
    if (!noPredictionFlag) spObjectPred <- sp::spTransform(x = spObjectPred, CRSobj = lonLatProjString)
  }

  ##################################
  # DEFINING CONTROL PARA.##########
  defaultControl <- list(Mlon = 1, Mlat = 1, Mtime = 1, randomSeed = 24, nuggetSD = 1e-5, numKnotsRes0 = 20L, J = 4L, numValuesForIS = 100, numIterOptim = 25L, distMethod = "haversine", normalHyperprior = TRUE, numISpropDistUpdates = 0, tipKnotsThinningRate = 1, credIntervalPercs = c(0.025, 0.975), timeJitterMaxInDays = 0, spaceJitterMax = 0)
  coordRanges <- .prepareCoordRanges(spObject = spObject, spObjectPred = spObjectPred, timePOSIXctVec = timePOSIXctVec, predTimePOSIXctVec = predTimePOSIXctVec)
  defaultControl <- c(defaultControl, coordRanges)

  for (i in names(control)) {
    defaultControl[[i]] <- control[[i]]
  }
  control <- defaultControl
  ##################################
  if (control$spaceJitterMax > 0) spObject@coords <- geoR::jitter2d(spObject@coords, max = control$spaceJitterMax)

  if (is.null(spatialRangeList)) {
    warning("Did not provide a starting value for the spatial range parameter. Using a fifth of the length of the training data bounding box, and setting hyperparameters mu = 'default starting value' and sigma = 'default starting value').")
    lowerLeftCorner <- c(control$lonRange[[1]], control$latRange[[1]])
    upperRightCorner <- c(control$lonRange[[2]], control$latRange[[2]])
    diagLengthInKm <- geosphere::distHaversine(p1 = lowerLeftCorner, p2 = upperRightCorner, r = 6378.137) # Distances are in kilometers.
    logRangePara <- log(diagLengthInKm/5)
    spatialRangeList <- list(start = logRangePara,
                             hyperpars = c(mu = logRangePara,
                                           sigma = logRangePara)
    )
  }

  if (is.null(timeRangeList)) {
    warning("Did not provide a starting value for the temporal range parameter. Using a fifth of the length of the time range (expressed in days), and setting hyperparameters mu = 'default starting value'  and sigma = 'default starting value'.")
    timeRangeList <- list(start = diff(control$timeRange)/5)
    timeRangeList$hyperpars <- c(mu = timeRangeList$start, sigma = timeRangeList$start)
  }

  if (!is.null(control$folderToSaveISpoints)) {
    numPoints <- length(list.files(path = control$folderToSaveISpoints, pattern = "ISoutput"))
    if (numPoints >= control$numValuesForIS) {
      control$IScompleted <- TRUE
    }
  }
  timeValues <- .ConvertPOSIXctInDays(timePOSIXctVec, min(c(timePOSIXctVec, predTimePOSIXctVec)))
  predTime <- NULL
  if (!noPredictionFlag) {
    predTime <- .ConvertPOSIXctInDays(predTimePOSIXctVec, min(timePOSIXctVec))
    predCovariateFrame <- as.matrix(predCovariateFrame)
    predCoords <- spObjectPred@coords
  } else {
    predCoords <- NULL
  }

  hyperStart <- .makeHyperStart(spatialRangeList = spatialRangeList, spatialSmoothnessList = spatialSmoothnessList, timeRangeList = timeRangeList, timeSmoothnessList = timeSmoothnessList, errorSDlist = errorSDlist, fixedEffSDlist = fixedEffSDlist, scaleList = scaleList)
  fixedHyperValues <- .makeFixedHyperValues(spatialRangeList = spatialRangeList, spatialSmoothnessList = spatialSmoothnessList, timeRangeList = timeRangeList, timeSmoothnessList = timeSmoothnessList, errorSDlist = errorSDlist, fixedEffSDlist = fixedEffSDlist, scaleList = scaleList)
  hyperpriorPars <- .makeHyperpriorPars(spatialRangeList = spatialRangeList, spatialSmoothnessList = spatialSmoothnessList, timeRangeList = timeRangeList, timeSmoothnessList = timeSmoothnessList, errorSDlist = errorSDlist, fixedEffSDlist = fixedEffSDlist, scaleList = scaleList)

  nestedGridsPointer <- setupNestedGrids(responseValues = responseVec, spCoords = spObject@coords, predCoords = predCoords, obsTime = timeValues, predTime = predTime, covariateMatrix = as.matrix(covariateFrame), predCovariateMatrix = predCovariateFrame, Mlon = control$Mlon, Mlat = control$Mlat, Mtime = control$Mtime, lonRange = control$lonRange, latRange = control$latRange, timeRange = control$timeRange, randomSeed = control$randomSeed, numKnotsRes0 = control$numKnotsRes0, J = control$J, distMethod = control$distMethod, MaternParsHyperpars = hyperpriorPars[c("space", "time", "scale")], fixedEffParsHyperpars = hyperpriorPars$fixedEffSD, errorParsHyperpars = hyperpriorPars$errorSD, FEmuVec = FEmuVec, nuggetSD = control$nuggetSD, normalHyperprior = control$normalHyperprior, tipKnotsThinningRate = control$tipKnotsThinningRate)$nestedGridsPointer

  # First we compute values relating to the hyperprior marginal distribution...

  computedValues <- .obtainISvalues(nestedGridsPointer = nestedGridsPointer, hyperStart = hyperStart, fixedHyperValues = fixedHyperValues, control = control)

  # Now, we obtain the marginal distribution of all mean parameters.
  cat("Computing moments for marginal posterior distributions...\n")
  computedValues$output <- .AddLogISweight(output = computedValues$output, distMode = computedValues$ISdistParas$mu, control = control)
  hyperMarginalMoments <- .ComputeHyperMarginalMoments(computedValues$output, control = control)
  FEmarginalMoments <- .ComputeFEmarginalMoments(computedValues$output, covNames = c("Intercept", colnames(covariateFrame)), control = control)

  outputList <- list(hyperMarginalMoments = hyperMarginalMoments$paraMoments, FEmarginalMoments = FEmarginalMoments, psiAndMargDistMatrix = hyperMarginalMoments$psiAndMargDistMatrix)

  if (!noPredictionFlag) {
    cat("Computing prediction moments... \n")
    outputList$predictionMoments <- .ComputeKrigingMoments(computedValues$output, nestedGridsPointer, control = control)
  }

  cat("Returning results... \n")
  outputList
}

.ConvertPOSIXctInDays <- function(POSIXctVec, baselinePOSIXct = 0) {
  as.numeric(POSIXctVec - baselinePOSIXct)/(3600 * 24)
}

.checkInputConsistency <- function(responseVec, covariateFrame, spatialCoordMat, timePOSIXctVec, predCovariateFrame, predSpatialCoordMat, predTimePOSIXctVec) {
  noPredictionFlag <- is.null(predCovariateFrame) | is.null(predSpatialCoordMat) | is.null(predTimePOSIXctVec)
  if (noPredictionFlag) {
    warning("Missing component for predictions: the model will be fitted, but no predictions will be produced.")
  }

  if (!identical(NULL, predCovariateFrame) | !identical(colnames(covariateFrame), colnames(predCovariateFrame))) {
    if (!identical(sort(colnames(covariateFrame)), sort(colnames(predCovariateFrame)))) {
      stop("Mismatch between covariates in training and test data. \n")
    }
  }
  if (!all(duplicated(c(length(responseVec), nrow(covariateFrame), nrow(spatialCoordMat), length(timePOSIXctVec)))[-1])) stop("Mismatch in the dimensions of responseVec, covariateFrame, spatialCoordMat, and/or timePOSIXctVec.")
  if (!noPredictionFlag) {
    if (!all(duplicated(c(nrow(predCovariateFrame), nrow(predSpatialCoordMat), length(predTimePOSIXctVec)))[-1])) stop("Mismatch in the dimensions of predCovariateFrame, predSpatialCoordMat, and/or predTimePOSIXctVec.")
  }
  NULL
}

.prepareCoordRanges <- function(spObject, spObjectPred, timePOSIXctVec, predTimePOSIXctVec) {
  spCoordRanges <- lapply(1:2, function(coordColIndex) {
    bufferSize <- 0.01
    coordinates <- spObject@coords[ , coordColIndex]
    combinedRangeNoBuffer <- range(coordinates)
    if (!is.null(spObjectPred)) {
      predCoordinates <- spObjectPred@coords[, coordColIndex]
      combinedRangeNoBuffer <- range(c(coordinates, predCoordinates))
    }
    combinedRangeNoBuffer + c(-bufferSize, bufferSize)
  })
  names(spCoordRanges) <- c("lonRange", "latRange")

  timeValuesRangeNoBuffer <- range(timePOSIXctVec)
  if (!is.null(spObjectPred)) timeValuesRangeNoBuffer <- range(c(timePOSIXctVec, predTimePOSIXctVec))
  timeBufferSize <- 10 # In seconds
  timeValuesRange <- timeValuesRangeNoBuffer + c(-timeBufferSize, timeBufferSize)
  timeRangeReshaped <- .ConvertPOSIXctInDays(timeValuesRange, min(timeValuesRangeNoBuffer))
  c(spCoordRanges, list(timeRange = timeRangeReshaped))
}

.makeFixedHyperValues <- function(spatialRangeList, spatialSmoothnessList, timeRangeList, timeSmoothnessList, errorSDlist, fixedEffSDlist, scaleList) {
  fixedHyper <- list()
  if (length(spatialRangeList) == 1) {
    fixedHyper$space <- c(rho = spatialRangeList[[1]])
  }
  if (length(spatialSmoothnessList) == 1) {
    fixedHyper$space <- c(fixedHyper$space["rho"], smoothness = spatialSmoothnessList[[1]])
    # 1e5 is used as a substitute for infinity, which is not understood by the C++ code.
    if (fixedHyper$space[["smoothness"]] > 1e5) fixedHyper$space[["smoothness"]] <- 1e5
  }
  if (length(timeRangeList) == 1) {
    fixedHyper$time <- c(rho = timeRangeList[[1]])
  }
  if (length(timeSmoothnessList) == 1) {
    fixedHyper$time <- c(fixedHyper$time["rho"], smoothness = timeSmoothnessList[[1]])
    if (fixedHyper$time[["smoothness"]] > 1e5) fixedHyper$time[["smoothness"]] <- 1e5
  }

  for (argName in c("errorSDlist", "fixedEffSDlist", "scaleList")) {
    if (length(get(argName)) == 1) {
      paraName <- substr(x = argName, start = 1, stop = nchar(argName) - 4)
      fixedHyper[[paraName]] <- get(argName)[[1]]
    }
  }
  fixedHyper
}

.makeHyperStart <- function(spatialRangeList, spatialSmoothnessList, timeRangeList, timeSmoothnessList, errorSDlist, fixedEffSDlist, scaleList) {
  hyperStart <- list()
  if (length(spatialRangeList) > 1) {
    hyperStart$space <- c(rho = spatialRangeList[[1]])
  }
  if (length(spatialSmoothnessList) > 1) {
    hyperStart$space <- c(hyperStart$space["rho"], smoothness = spatialSmoothnessList[[1]])
    if (hyperStart$space[["smoothness"]] > 1e5) hyperStart$space[["smoothness"]] <- 1e5
  }
  if (length(timeRangeList) > 1) {
    hyperStart$time <- c(rho = timeRangeList[[1]])
  }
  if (length(timeSmoothnessList) > 1) {
    hyperStart$time <- c(hyperStart$time["rho"], smoothness = timeSmoothnessList[[1]])
    if (hyperStart$time[["smoothness"]] > 1e5) hyperStart$time[["smoothness"]] <- 1e5
  }

  for (argName in c("errorSDlist", "fixedEffSDlist", "scaleList")) {
    if (length(get(argName)) > 1) {
      paraName <- substr(x = argName, start = 1, stop = nchar(argName) - 4)
      hyperStart[[paraName]] <- get(argName)[[1]]
    }
  }
  hyperStart
}

.makeHyperpriorPars <- function(spatialRangeList, spatialSmoothnessList, timeRangeList, timeSmoothnessList, errorSDlist, fixedEffSDlist, scaleList) {
  hyperpriorPars <- list(
    space = list(rho = c(mu = 0, sigma = 0), smoothness = c(mu = 0, sigma = 0)),
    time = list(rho = c(mu = 0, sigma = 0), smoothness = c(mu = 0, sigma = 0)),
    errorSD = c(mu = 0, sigma = 0),
    fixedEffSD = c(mu = 0, sigma = 0),
    scale = c(mu = 0, sigma = 0)
  )
  if (length(spatialRangeList) > 1) {
    hyperpriorPars$space$rho <- spatialRangeList[[2]]
  }
  if (length(spatialSmoothnessList) > 1) {
    hyperpriorPars$space$smoothness <- spatialSmoothnessList[[2]]
  }
  if (length(timeRangeList) > 1) {
    hyperpriorPars$time$rho <- timeRangeList[[2]]
  }
  if (length(timeSmoothnessList) > 1) {
    hyperpriorPars$time$smoothness <- timeSmoothnessList[[2]]
  }

  for (argName in c("errorSDlist", "fixedEffSDlist", "scaleList")) {
    if (length(get(argName)) > 1) {
      paraName <- substr(x = argName, start = 1, stop = nchar(argName) - 4)
      hyperpriorPars[[paraName]] <- get(argName)[[2]]
    }
  }
  hyperpriorPars
}

.AddLogISweight <- function(output, distMode, control) {
  adaptiveISphaseVector <- rep(1:(control$numISpropDistUpdates + 1), each = ceiling(control$numValuesForIS / (control$numISpropDistUpdates + 1)))[1:min(control$numValuesForIS, length(output))]
  for (phase in 1:max(adaptiveISphaseVector)) {
    itersInPhase <- which(adaptiveISphaseVector == phase)
    weightModifs <- sapply(output[itersInPhase], FUN = function(outputElement) {
      mvtnorm::dmvnorm(x = outputElement$x, mean = distMode, sigma = outputElement$varCovar, log = TRUE)
    })
    discreteLogJointValues <- sapply(output[itersInPhase], '[[', "logJointValue")
    logWeights <- discreteLogJointValues - weightModifs - log(length(discreteLogJointValues))
    maxLogWeights <- max(logWeights)
    logPropConstantIS <- maxLogWeights + log(sum(exp(logWeights - maxLogWeights)))
    logStandardisedWeights <- logWeights - logPropConstantIS
    for (j in seq_along(itersInPhase)) {
      output[[itersInPhase[[j]]]]$logISweight <- logStandardisedWeights[[j]]
    }
  }
  output
}

.obtainISvalues <- function(nestedGridsPointer, hyperStart, fixedHyperValues, control) {
  iterCounter <- 0
  funForOptim <- function(xOnLogScale, namesXstartValues) {
    iterCounter <<- iterCounter + 1
    cat("Performing evaluation ", iterCounter, ".\n")
    names(xOnLogScale) <- namesXstartValues
    xTrans <- exp(xOnLogScale)
    unlistedFixedHyperValues <- unlist(fixedHyperValues)
    if (control$normalHyperprior) {
      unlistedFixedHyperValues <- exp(unlistedFixedHyperValues)
    }
    hyperList <- .prepareHyperList(xTrans, fixedHyperValuesUnlisted = unlistedFixedHyperValues)
    returnedValue <- -.LogJointHyperMarginal(treePointer = nestedGridsPointer, hyperparaValues = hyperList, recordFullConditional = FALSE, processPredictions = FALSE)
    returnX <- xOnLogScale
    if (!control$normalHyperprior) returnX <- exp(xOnLogScale)
    returnedValue
  }

  gradForOptim <- function(xOnLogScale, namesXstartValues) {
    names(xOnLogScale) <- namesXstartValues
    numDeriv::grad(func = funForOptim, x = xOnLogScale, method = "simple", namesXstartValues = namesXstartValues)
  }
  # A short optimisation first...

  # Remember that if the priors for hyperparameters are normal, starting values are given on the log-scale.
  xStartValues <- unlist(hyperStart)

  lowerBound <- rep(1e-10, length(xStartValues))
  if (control$normalHyperprior) {
    lowerBound <- log(lowerBound)
  }
  upperBound <- rep(Inf, length(xStartValues))
  if (!is.null(control$upperBound)) {
    upperBound <- unlist(control$upperBound)
  }
  names(upperBound) <- names(lowerBound) <- names(xStartValues)
  cat("Optimising... \n")
  if (!tryCatch(file.exists(control$fileToSaveOptOutput), error = function(e) FALSE)) { # The tryCatch is necessary to ensure that an error does not occur if control$fileToSaveOptOutput is NULL. If it is undefined, we want the optimisation to take place.
    x0val <- xStartValues
    if (!control$normalHyperprior) {
      x0val <- log(xStartValues)
    }
    opt <- nloptr::lbfgs(x0 = x0val, lower = lowerBound, upper = upperBound, fn = funForOptim, gr = gradForOptim, control = list(xtol_rel = 1e-3, maxeval = control$numIterOptim), namesXstartValues = names(xStartValues))
    solution <- opt$par
    if (!control$normalHyperprior) {
      solution <- exp(opt$par)
    }
    cat("Computing precision matrix at mode... \n")
    precisionMat <- numDeriv::hessian(func = funForOptim, x = solution, namesXstartValues = names(xStartValues))
    varCovar <- tryCatch(expr = solve(precisionMat), error = function(e) e)

    if ("error" %in% class(varCovar)) {
      warning("Hessian is singular! Adding nugget along diagonal to correct...", immediate. = TRUE)
      varCovar <- solve(precisionMat + 1e-10 * diag(nrow(precisionMat)))
    }
    eigenVarCovar <- eigen(varCovar)
    if (any(eigenVarCovar$values < 0)) {
      warning("Covariance matrix for proposal distribution is not positive definite! Correcting with Matrix::nearPD.", immediate. = TRUE)
      varCovar <- as.matrix(Matrix::nearPD(x = varCovar, ensureSymmetry = TRUE)$mat)
    }

    if (!is.null(control$fileToSaveOptOutput)) {
      tryCatch(expr = save(opt, file = control$fileToSaveOptOutput), error = function(e) {warning("Could not save optimisation results! It will not be possible to resume the function after the optimisation step, should you need to interrupt the current run. Check the file name you provided.", immediate. = TRUE)})
      filenameForVarCovar <- paste(substr(control$fileToSaveOptOutput, start = 1, stop = gregexpr(pattern = ".Rdata", text = control$fileToSaveOptOutput)[[1]] - 1), "_ISvarCovar.Rdata", sep = "")
      tryCatch(expr =  save(varCovar, file = filenameForVarCovar), error = function(e) invisible(NULL))
    }
    cat("Optimised values:", solution)
  } else {
    load(control$fileToSaveOptOutput)
    solution <- opt$par
    if (!control$normalHyperprior) {
      solution <- exp(opt$par)
    }
    varCovarFilename <- paste(substr(control$fileToSaveOptOutput, start = 1, stop = gregexpr(pattern = ".Rdata", text = control$fileToSaveOptOutput)[[1]] - 1), "_ISvarCovar.Rdata", sep = "")
    loadName <- load(varCovarFilename) # Restores covariance matrix
    varCovar <- get(loadName)
    rm(loadName)
  }

  if (!is.null(control$envirForTest)) {
    assign(x = "Hmat", value = GetHmat(nestedGridsPointer), envir = control$envirForTest)
  }

  opt$value <- -opt$value # Correcting for the inversion used to maximise instead of minimise

  cat("Running IS algorithm... \n")
  ISvaluesList <- .ISfct(distMode = solution, ISvarCovar = varCovar, nestedGridsPointer = nestedGridsPointer, namesXstartValues = names(xStartValues), fixedHyperValues = fixedHyperValues,  control = control)
  cat("IS algorithm completed... \n")

  keepIndices <- sapply(ISvaluesList$output, function(x) class(x$logJointValue) == "numeric")
  ISvaluesList$output <- ISvaluesList$output[keepIndices]
  ISvaluesList$ISdistParas <- list(mu = solution, cov = varCovar)
  ISvaluesList
}

.prepareHyperList <- function(hyperStartUnlisted, fixedHyperValuesUnlisted) {
  paraValues <- sapply(c("error", "fixed", "space.rho", "space.smoothness", "time.rho", "time.smoothness", "scale"), function(paraName) {
    if (length(pos <- grep(pattern = paraName, x = names(hyperStartUnlisted))) > 0) {
      argValue <- hyperStartUnlisted[[pos]]
    } else if (length(pos <- grep(pattern = paraName, x = names(fixedHyperValuesUnlisted)))) {
      argValue <- fixedHyperValuesUnlisted[[pos]] # Already on the correct scale, ie not log scale
    } else {
      stop(paste("Missing hyperparameter specification: ", paraName, "! \n", sep = ""))
    }
    argValue
  })
  list(space = c(rho = paraValues[["space.rho"]], smoothness = paraValues[["space.smoothness"]]), time = c(rho = paraValues[["time.rho"]], smoothness = paraValues[["time.smoothness"]]), scale = paraValues[["scale"]], errorSD = paraValues[["error"]], fixedEffSD = paraValues[["fixed"]])
}

.ISfct <- function(distMode, ISvarCovar, nestedGridsPointer, namesXstartValues, fixedHyperValues, control) {
  updatedISvarCovar <- ISvarCovar # Cov. matrix which was updated with adaptive IS.
  baseNumItersInPhase <- ceiling(control$numValuesForIS/(control$numISpropDistUpdates + 1))
  numPhases <- control$numISpropDistUpdates + 1
  numItersInPhaseVec <- rep(baseNumItersInPhase, numPhases)
  numItersInPhaseVec[length(numItersInPhaseVec)] <- control$numValuesForIS - baseNumItersInPhase * (numPhases - 1)
  generalCounter <- 0
  output <- vector("list", length = control$numValuesForIS)
  startingPhase <- currentIterInPhase <- 1
  if (tryCatch(dir.exists(control$folderToSaveISpoints), error = function(e) FALSE)) {
    cat("Loading previously processed IS points... \n")
    loadResult <- lapply(1:numPhases, function(phaseNumber) {
      lapply(list.files(control$folderToSaveISpoints, full.names = TRUE, pattern = paste("ISoutputPhase", phaseNumber, sep = "")), function(x) get(load(x)))
    })

    phaseIndices <- sapply(1:numPhases, function(phaseNumber) {
      filesToLoad <- list.files(control$folderToSaveISpoints, full.names = TRUE, pattern = paste("ISoutputPhase", phaseNumber, sep = ""))
      length(filesToLoad)
    })

    if (numPhases > 1) {
      startingPhase <- max(match(0, phaseIndices) - 1, 1)
    }
    currentIterInPhase <- phaseIndices[[startingPhase]] + 1

    numLoadedResults <- sum(sapply(loadResult, FUN = length))

    if (numLoadedResults > 0) {
      output[1:numLoadedResults] <- unlist(loadResult, recursive = FALSE)
    }
    if (identical(control$IScompleted, TRUE) | (numLoadedResults == control$numValuesForIS)) { # If we indicate that the IS is completed, the function will simply stop sampling points and use saved points only.
      output <- output[1:numLoadedResults]
      return(list(output = output))
    }
    generalCounter <- numLoadedResults
  }

  for (phase in startingPhase:numPhases) {
    numItersInPhase <- numItersInPhaseVec[[phase]]

    if (!control$normalHyperprior) {
      smallMVR <- function() {
        container <- NULL
        repeat {
          container <- drop(mvtnorm::rmvnorm(n = 1, mean = distMode, sigma = ISvarCovar))
          if (all(container > 0)) break
        }
        container
      }
      paraGrid <- t(replicate(n = numItersInPhase - currentIterInPhase + 1, expr = smallMVR()))
    } else {
      paraGrid <- mvtnorm::rmvnorm(n = numItersInPhase - currentIterInPhase + 1, mean = distMode, sigma = updatedISvarCovar)
    }
    colnames(paraGrid) <- namesXstartValues

    for (i in seq_along(currentIterInPhase:numItersInPhase)) {
      generalCounter <- generalCounter + 1
      cat("Processing grid value ", generalCounter, "... \n")
      xVec <- paraGrid[i, ]
      names(xVec) <- colnames(paraGrid)
      if (control$normalHyperprior) {
        xVec <- exp(xVec)
      }
      output[[generalCounter]] <- .funForGridEst(xNonLogScale = xVec, treePointer = nestedGridsPointer, fixedHyperValues = fixedHyperValues, computePrediction = TRUE, control = control)
      output[[generalCounter]]$varCovar <- updatedISvarCovar
      if (!is.null(control$folderToSaveISpoints)) {
        if (!dir.exists(control$folderToSaveISpoints)) {
          tryCatch(expr = dir.create(path = control$folderToSaveISpoints), error = function(e) warning("Could not create directory to save importance sampling results! Check the directory name you provided. It will not be possible to resume the IS algorithm should you decide to interrupt this run.", immediate. = TRUE))
        }
        filename <- paste(control$folderToSaveISpoints, "/", "ISoutputPhase", phase,"_iter", generalCounter, ".Rdata", sep = "")
        objectToSave <- output[[generalCounter]]
        tryCatch(expr = save(objectToSave, file = filename, compress = TRUE), error = function(e) invisible(NULL))

        if (generalCounter == 1) { # On the first iteration, the vector to restore the order of predictions must be saved to allow for the predicted values to be re-ordered after a resume in which CreateHmatrixPred is not called. This situation occurs when we specify control$IScompleted=TRUE.
          predOrder <- GetPredObsOrder(nestedGridsPointer)
          filenameForPredOrder <- paste(control$folderToSaveISpoint, "/predObsOrder.Rdata", sep = "")
          tryCatch(expr = save(predOrder, file = filenameForPredOrder, compress = TRUE), error = function(e) invisible(NULL))
        }
      }
    }
    currentIterInPhase <- 1 # currentIterInPhase is there to allow the function to resume after interruption
    if (control$numISpropDistUpdates > 0) {
      outputIndices <- ((phase - 1) * numItersInPhase + 1):(min(control$numValuesForIS, phase * numItersInPhase))
      weightModifs <- sapply(output[outputIndices], MARGIN = 1, FUN = function(outputElement) {
        mvtnorm::dmvnorm(x = outputElement$x, mean = distMode, sigma = updatedISvarCovar, log = TRUE)
      })
      discreteLogJointValues <- sapply(output[outputIndices], '[[', "logJointValue")
      logWeights <- discreteLogJointValues - weightModifs - log(length(discreteLogJointValues))
      maxLogWeights <- max(logWeights)
      logPropConstantIS <- maxLogWeights + log(sum(exp(logWeights - maxLogWeights)))
      logStandardisedWeights <- logWeights - logPropConstantIS
      updatedISvarCovar <- cov.wt(x = paraGrid, wt = exp(logStandardisedWeights), center = distMode)$cov

      if (any(eigen(updatedISvarCovar)$values < 0)) {
        updatedISvarCovar <- as.matrix(Matrix::nearPD(x = updatedISvarCovar, ensureSymmetry = TRUE)$mat)
      }
      if (Inf %in% updatedISvarCovar) {
        stop("Only one weight for the estimation of the covariance matrix. Stop here for now.")
      }
      if (!is.null(control$fileToSaveOptOutput)) {
        varCovarFilename <- paste(substr(control$fileToSaveOptOutput, start = 1, stop = gregexpr(pattern = ".Rdata", text = control$fileToSaveOptOutput)[[1]] - 1), "_ISvarCovar.Rdata", sep = "")
        tryCatch(expr = save(updatedISvarCovar, file = varCovarFilename), error = function(e) warning("Could not save proposal variance-covariance matrix in IS run! Check control$fileToSaveOptOutput.", immediate. = TRUE)) # We update the matrix in memory. If the function is restarted, it will resume with the last matrix saved, which is what we want.
      }
    }
  }

  list(output = output)
}

.funForGridEst <- function(xNonLogScale, treePointer, fixedHyperValues, computePrediction, control) {
  fixedHyperValuesUnlisted <- unlist(fixedHyperValues)
  if (control$normalHyperprior) {
    fixedHyperValuesUnlisted <- exp(fixedHyperValuesUnlisted)
  }
  hyperList <- .prepareHyperList(hyperStartUnlisted = xNonLogScale, fixedHyperValuesUnlisted = fixedHyperValuesUnlisted)

  # logJointValue <- tryCatch(expr = .LogJointHyperMarginal(treePointer = treePointer, hyperparaValues = hyperList, recordFullConditional = FALSE, processPredictions = TRUE), error = function(e) e)
  cat("Processed value: \n")
  print(hyperList)
  logJointValue <- .LogJointHyperMarginal(treePointer = treePointer, hyperparaValues = hyperList, recordFullConditional = FALSE, processPredictions = TRUE)
  x <- xNonLogScale
  if (control$normalHyperprior) {
    x <- log(xNonLogScale)
    hyperList$fixedEffSD <- log(hyperList$fixedEffSD)
    hyperList$errorSD <- log(hyperList$errorSD)
    hyperList$space <- log(hyperList$space)
    hyperList$time <- log(hyperList$time)
    hyperList$scale <- log(hyperList$scale)
  }
  aList <- list(x = x, fixedEffSD = hyperList$fixedEffSD , errorSD = hyperList$errorSD, MaternHyperpars = hyperList[c("space", "time", "scale")], logJointValue = logJointValue)
  if (computePrediction) {
    aList$CondPredStats <- ComputeCondPredStats(treePointer)
  }
  # Running .LogJointHyperMarginal stores in the tree pointed by nestedGridsPointer the full conditional mean and SDs when recordFullConditional = TRUE. We can get them with the simple functions I call now.
  aList$FullCondMean <- GetFullCondMean(treePointer)
  aList$FullCondSDs <- GetFullCondSDs(treePointer)
  aList
}

.ComputeHyperMarginalMoments <- function(hyperparaList, control) {
  domainCheck <- sapply(hyperparaList, function(x) x$logJointValue > -Inf)
  hyperparaList <- hyperparaList[domainCheck]
  psiAndMargDistMatrix <- t(sapply(seq_along(hyperparaList), function(hyperparaIndex) c(unlist(hyperparaList[[hyperparaIndex]]$MaternHyperpars), fixedEffSD = hyperparaList[[hyperparaIndex]]$fixedEffSD, errorSD = hyperparaList[[hyperparaIndex]]$errorSD, logJointValue = hyperparaList[[hyperparaIndex]]$logJointValue, ISweight = exp(hyperparaList[[hyperparaIndex]]$logISweight))))
  rownames(psiAndMargDistMatrix) <- NULL
  adaptiveISphaseVector <- rep(1:(control$numISpropDistUpdates + 1), each = ceiling(control$numValuesForIS / (control$numISpropDistUpdates + 1)))[1:min(control$numValuesForIS, length(hyperparaList))]
  funToGetParaMoments <- function(hyperparaIndex) {
    # meanValue <- sum(psiAndMargDistMatrix[, hyperparaIndex] * psiAndMargDistMatrix[, ncol(psiAndMargDistMatrix)])
    # sdValue <- sqrt(sum(psiAndMargDistMatrix[, hyperparaIndex]^2 * psiAndMargDistMatrix[, ncol(psiAndMargDistMatrix)]) - meanValue^2)
    meanValue <- .adaptiveIS(x = psiAndMargDistMatrix[, hyperparaIndex], ISweights = psiAndMargDistMatrix[, "ISweight"], phaseVector = adaptiveISphaseVector)

    sdValue <- sqrt(.adaptiveIS(x = psiAndMargDistMatrix[, hyperparaIndex]^2, ISweights = psiAndMargDistMatrix[, "ISweight"], phaseVector = adaptiveISphaseVector) - meanValue^2)

    skewnessValue <- .adaptiveIS(x = (psiAndMargDistMatrix[, hyperparaIndex] - meanValue)^3/sdValue^2, ISweights = psiAndMargDistMatrix[, "ISweight"], phaseVector = adaptiveISphaseVector)
    credIntBounds <- list(bounds = c(NA, NA), alpha = NA, omega = NA , xi = NA, delta = NA)
    if (!((sdValue == 0) | is.na(sdValue))) {
      credIntBounds <- .ComputeCredIntervalSkewNorm(control$credIntervalPercs, meanValue = meanValue, sdValue = sdValue, skewnessValue = skewnessValue)
    }
    c(mean = meanValue, StdDev = sdValue, skewness = skewnessValue, credIntBounds$bounds)
  }
  paraMoments <- t(sapply(1:(ncol(psiAndMargDistMatrix) - 2), FUN = funToGetParaMoments))
  CInames <- paste("CredInt_", round(control$credIntervalPercs, 3)*100, "%", sep = "")
  colnames(paraMoments) <- c("Mean", "StdDev", "Skewness", CInames)
  rownames(paraMoments) <- head(colnames(psiAndMargDistMatrix), n = -2)
  list(paraMoments = as.data.frame(paraMoments), psiAndMargDistMatrix = psiAndMargDistMatrix)
}

.ComputeCredIntervalSkewNorm <- function(p = c(0.025, 0.975), meanValue, sdValue, skewnessValue) {
  skewNormalDelta <- sign(skewnessValue) * sqrt(
    pi/2 * abs(skewnessValue)^(2/3) /
      (abs(skewnessValue)^(2/3) + ((4 - pi)/2)^(2/3))
  )
  skewNormalAlpha <- skewNormalDelta/sqrt(1 - skewNormalDelta^2)
  skewNormalOmega <- sdValue / sqrt(1 - 2 * skewNormalDelta^2/pi)
  skewNormalXi <- meanValue - skewNormalOmega * skewNormalDelta * sqrt(2/pi)
  bounds <- sn::qsn(p = p, xi = skewNormalXi, omega = skewNormalOmega, alpha = skewNormalAlpha)
  names(bounds) <- paste("CredInt_", round(p, 3)*100, "%", sep = "")
  list(bounds = bounds, alpha = skewNormalAlpha, omega = skewNormalOmega, xi = skewNormalXi, delta = skewNormalDelta)
}

.ComputeFEmarginalMoments <- function(hyperparaList, covNames, control) {
  adaptiveISphaseVector <- rep(1:(control$numISpropDistUpdates + 1), each = ceiling(control$numValuesForIS / (control$numISpropDistUpdates + 1)))[1:min(control$numValuesForIS, length(hyperparaList))]
  logISweightVector <- sapply(hyperparaList, function(x) x$logISweight)
  marginalMeans <- sapply(1:length(covNames), function(paraIndex) {
    meanVector <- sapply(hyperparaList, function(x) x$FullCondMean[[paraIndex]])
    .adaptiveIS(x = meanVector, ISweights = exp(logISweightVector), phaseVector = adaptiveISphaseVector)
  })
  marginalSecondMoments <- sapply(1:length(covNames), function(paraIndex) {
    meanVector <- sapply(hyperparaList, function(x) x$FullCondMean[[paraIndex]])
    sdVector <- sapply(hyperparaList, function(x) x$FullCondSDs[[paraIndex]])
    secondMomentVec <- sdVector^2 + meanVector^2
    .adaptiveIS(x = secondMomentVec, ISweights = exp(logISweightVector), phaseVector = adaptiveISphaseVector)
  })
  marginalSDs <- sqrt(marginalSecondMoments - marginalMeans^2)
  credIntFrame <- ComputeFEcredInts(control$credIntervalPercs, hyperparaList, marginalMeans, marginalSDs)
  outputFrame <- cbind(data.frame(Mean = marginalMeans, StdDev = marginalSDs), credIntFrame)
  rownames(outputFrame) <- covNames
  outputFrame
}

ComputeFEcredInts <- function(p = c(0.025, 0.975), hyperparaList, marginalMeans, marginalSDs) {
  valuesRanges <- lapply(1:length(marginalMeans), function(x) {
    c(-3 * marginalSDs[[x]], 3 * marginalSDs[[x]]) + marginalMeans[[x]]
  })
  distValuesByFEpar <- lapply(seq_along(marginalMeans), function(FEindex) {
    valuesToConsider <- seq(
      from = valuesRanges[[FEindex]][1],
      to = valuesRanges[[FEindex]][2],
      length.out = 200
    )
    distValuesByISiter <- sapply(hyperparaList, function(hyperparaListElement) {
      dnorm(valuesToConsider, mean = hyperparaListElement$FullCondMean[[FEindex]], sd = hyperparaListElement$FullCondSDs[[FEindex]]) * exp(hyperparaListElement$logISweight)
    })
    summedValues <- Reduce("+", distValuesByISiter)
    data.frame(x = valuesToConsider, values = rowSums(distValuesByISiter))
  })
  boundsByFEpar <- lapply(distValuesByFEpar, FUN = function(distFrame) {
    distributionValuesNormalised <- distFrame$values/sum(distFrame$values)
    DFvalues <- cumsum(distributionValuesNormalised)
    leftBoundPos <- match(TRUE, DFvalues >= p[1])
    rightBoundPos <- match(TRUE, DFvalues >= p[2])
    c(distFrame$x[[leftBoundPos]], distFrame$x[[rightBoundPos]])
  })
  boundsFrame <- as.data.frame(do.call("rbind", boundsByFEpar))
  colnames(boundsFrame) <- paste("CredInt_", round(p, 3)*100, "%", sep = "")
  boundsFrame
}

.ComputeKrigingMoments <- function(hyperparaList, treePointer, control) {
  adaptiveISphaseVector <- rep(1:(control$numISpropDistUpdates + 1), each = ceiling(control$numValuesForIS / (control$numISpropDistUpdates + 1)))[1:min(control$numValuesForIS, length(hyperparaList))]
  logISweightVector <- sapply(hyperparaList, function(x) x$logISweight)
  krigingMeans <- sapply(1:length(hyperparaList[[1]]$CondPredStats$Hmean), function(predObsIndex) {
    predVector <- sapply(hyperparaList, function(x) x$CondPredStats$Hmean[[predObsIndex]])
    .adaptiveIS(x = predVector, ISweights = exp(logISweightVector), phaseVector = adaptiveISphaseVector)
  })

  varE <- sapply(1:length(hyperparaList[[1]]$CondPredStats$Hmean), function(predObsIndex) {
    predVector <- sapply(hyperparaList, function(x) x$CondPredStats$Hmean[[predObsIndex]]^2)
    .adaptiveIS(x = predVector, ISweights = exp(logISweightVector), phaseVector = adaptiveISphaseVector) - krigingMeans[[predObsIndex]]^2
  })

  Evar <- sapply(1:length(hyperparaList[[1]]$CondPredStats$Hmean), function(predObsIndex) {
    EvarVector <- sapply(hyperparaList, function(x) x$CondPredStats$Evar[[predObsIndex]]^2)
    .adaptiveIS(x = EvarVector, ISweights = exp(logISweightVector), phaseVector = adaptiveISphaseVector)
  })

  predObsOrder <- GetPredObsOrder(treePointer = treePointer)
  if (identical(control$IScompleted, TRUE)) {
    predObsOrder <- tryCatch(expr = get(load(paste(control$folderToSaveISpoints, "/predObsOrder.Rdata", sep = ""))), error = function(e) stop("Could not load prediction order file (named predObsOrder.Rdata). It is not in control$folderToSaveISpoints. Cannot produce summary results... Exiting."))
  }
  data.frame(predictMeans = krigingMeans[order(predObsOrder)], predictSDs = sqrt(varE + Evar)[order(predObsOrder)])
}

.LogJointHyperMarginal <- function(treePointer, hyperparaValues, recordFullConditional, processPredictions = FALSE) {
  LogJointHyperMarginalToWrap(treePointer = treePointer, MaternHyperpars = hyperparaValues[c("space", "time", "scale")], fixedEffSD = hyperparaValues$fixedEffSD, errorSD = hyperparaValues$errorSD, recordFullConditional = TRUE, processPredictions = processPredictions)
}

.ComputeLogJointHyperMarginal <- function(hyperparaMatrix, spacetimeData, predictionData, hyperStart, fixedHyperValues, hyperpriorPars, FEmuVec, control) {
  # CHECKS #########################
  if (!is.null(predictionData)) {
    if (!identical(colnames(spacetimeData@data)[-1], colnames(predictionData@data))) {
      stop("Mismatch between covariates in training and test data. \n")
    }
  }
  ##################################
  # DEFINING CONTROL PARA.##########
  defaultControl <- list(Mlon = 1, Mlat = 1, Mtime = 1, randomSeed = 24, nuggetSD = 0.00001, numKnotsRes0 = 20L, J = 4L, numValuesForIS = 200, numIterOptim = 200L, distMethod = "haversine", normalHyperprior = FALSE)
  if (length(position <- grep(colnames(spacetimeData@sp@coords), pattern = "lon")) >= 1) {
    colnames(spacetimeData@sp@coords)[[position[[1]]]] <- "x"
    if (!is.null(predictionData)) {
      colnames(predictionData@sp@coords)[[position[[1]]]] <- "x"
    }
  }
  if (length(position <- grep(colnames(spacetimeData@sp@coords), pattern = "lat")) >= 1) {
    colnames(spacetimeData@sp@coords)[[position[[1]]]] <- "y"
    if (!is.null(predictionData)) {
      colnames(predictionData@sp@coords)[[position[[1]]]] <- "y"
    }
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
  ##################################
  # 1e5 is used as a substitute for infinity, which is not understood by the C++ code.
  if ("smoothness" %in% names(hyperStart$space)) {
    if (hyperStart$space[["smoothness"]] > 1e5) hyperStart$space[["smoothness"]] <- 1e5
  } else {
    if (fixedHyperValues$space[["smoothness"]] > 1e5) fixedHyperValues$space[["smoothness"]] <- 1e5
  }

  if ("smoothness" %in% names(hyperStart$time)) {
    if (hyperStart$time[["smoothness"]] > 1e5) hyperStart$time[["smoothness"]] <- 1e5
  } else {
    if (fixedHyperValues$time[["smoothness"]] > 1e5) fixedHyperValues$time[["smoothness"]] <- 1e5
  }

  for (i in names(control)) {
    defaultControl[[i]] <- control[[i]]
  }

  control <- defaultControl

  if (!is.null(control$folderToSaveISpoints)) {
    numPoints <- length(list.files(path = control$folderToSaveISpoints, pattern = "ISoutput"))
    if (numPoints >= control$numValuesForIS) {
      control$IScompleted <- TRUE
    }
  }
  dataCoordinates <- spacetimeData@sp@coords
  predCoordinates <- predictionData@sp@coords
  predCovariates <- as.matrix(predictionData@data)

  timeRangeReshaped <- as.integer(control$timeRange)/(3600*24)
  timeBaseline <- min(timeRangeReshaped)
  timeValues <- as.integer(time(spacetimeData@time))/(3600*24) - timeBaseline # The division is to obtain values in days.
  predTime <- as.integer(time(predictionData))/(3600*24) - timeBaseline
  timeRangeReshaped <- timeRangeReshaped - timeBaseline

  covariateMatrix <- as.matrix(spacetimeData@data[, -1, drop = FALSE])

  nestedGridsPointer <- setupNestedGrids(responseValues = spacetimeData@data[, 1], spCoords = dataCoordinates, predCoords = predCoordinates, obsTime = timeValues, predTime = predTime, covariateMatrix = covariateMatrix, predCovariateMatrix = predCovariates, Mlon = control$Mlon, Mlat = control$Mlat, Mtime = control$Mtime, lonRange = control$lonRange, latRange = control$latRange, timeRange = timeRangeReshaped, randomSeed = control$randomSeed, numKnotsRes0 = control$numKnotsRes0, J = control$J, distMethod = control$distMethod, MaternParsHyperpars = hyperpriorPars[c("space", "time", "scale")], fixedEffParsHyperpars = hyperpriorPars$fixedEffSD, errorParsHyperpars = hyperpriorPars$errorSD, FEmuVec = FEmuVec, nuggetSD = control$nuggetSD, normalHyperprior = control$normalHyperprior, tipKnotsThinningRate = control$tipKnotsThinningRate)$nestedGridsPointer

  funForOptim <- function(xOnLogScale, namesXstartValues) {
    names(xOnLogScale) <- namesXstartValues
    xTrans <- exp(xOnLogScale)
    unlistedFixedHyperValues <- unlist(fixedHyperValues)
    if (control$normalHyperprior) {
      unlistedFixedHyperValues <- exp(unlistedFixedHyperValues)
    }
    hyperList <- MRAinla:::.prepareHyperList(xTrans, fixedHyperValuesUnlisted = unlistedFixedHyperValues)
    returnedValue <- -.LogJointHyperMarginal(treePointer = nestedGridsPointer, hyperparaValues = hyperList, recordFullConditional = FALSE, processPredictions = FALSE)
    returnedValue
  }

  gradForOptim <- function(xOnLogScale, namesXstartValues) {
    names(xOnLogScale) <- namesXstartValues
    numDeriv::grad(func = funForOptim, x = xOnLogScale, method = "simple", namesXstartValues = namesXstartValues)
  }
  # A short optimisation first...

  # Remember that if the priors for hyperparameters are normal, starting values are given on the log-scale.
  xStartValues <- unlist(hyperStart)

  lowerBound <- rep(1e-10, length(xStartValues))
  if (control$normalHyperprior) {
    lowerBound <- log(lowerBound)
  }
  upperBound <- rep(Inf, length(xStartValues))
  if (!is.null(control$upperBound)) {
    upperBound <- unlist(control$upperBound)
  }
  names(upperBound) <- names(lowerBound) <- names(xStartValues)
  cat("Optimising... \n")
  if (!tryCatch(file.exists(control$fileToSaveOptOutput), error = function(e) FALSE)) { # The tryCatch is necessary to ensure that an error does not occur if control$fileToSaveOptOutput is NULL. If it is undefined, we want the optimisation to take place.
    x0val <- xStartValues
    if (!control$normalHyperprior) {
      x0val <- log(xStartValues)
    }
    opt <- nloptr::lbfgs(x0 = x0val, lower = lowerBound, upper = upperBound, fn = funForOptim, gr = gradForOptim, control = list(xtol_rel = 1e-3, maxeval = control$numIterOptim), namesXstartValues = names(xStartValues))
    solution <- opt$par
    if (!control$normalHyperprior) {
      solution <- exp(opt$par)
    }
  } else {
    load(control$fileToSaveOptOutput)
    solution <- opt$par
    if (!control$normalHyperprior) {
      solution <- exp(opt$par)
    }
    varCovarFilename <- paste(substr(control$fileToSaveOptOutput, start = 1, stop = gregexpr(pattern = ".Rdata", text = control$fileToSaveOptOutput)[[1]] - 1), "_ISvarCovar.Rdata", sep = "")
    load(varCovarFilename) # Restores varCovar
  }

  hyperparasFormatted <- lapply(1:nrow(hyperparaMatrix), function(lineIndex) {
    xTrans <- hyperparaMatrix[lineIndex, ]
    unlistedFixedHyperValues <- unlist(fixedHyperValues)
    if (control$normalHyperprior) {
      xTrans <- exp(hyperparaMatrix[lineIndex, ])
      unlistedFixedHyperValues <- exp(unlistedFixedHyperValues)
    }
    hyperList <- MRAinla:::.prepareHyperList(xTrans, fixedHyperValuesUnlisted = unlistedFixedHyperValues)
  })

  sapply(hyperparasFormatted, FUN = .LogJointHyperMarginal, treePointer = nestedGridsPointer, recordFullConditional = FALSE, processPredictions = FALSE)
}

# See http://statweb.stanford.edu/~owen/pubtalks/AdaptiveISweb.pdf

.adaptiveIS <- function(x, ISweights, phaseVector) {
  resultPerPhase <- sapply(1:max(phaseVector), function(phaseNum) {
    itersToConsider <- phaseVector == phaseNum
    sum(x[itersToConsider] * ISweights[itersToConsider])
  })
  adaptiveISweights <- sqrt(1:max(phaseVector))
  sum(resultPerPhase * adaptiveISweights)/sum(adaptiveISweights)
}
