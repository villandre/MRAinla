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

MRA_INLA <- function(spacetimeData, errorSDstart, fixedEffSDstart, MRAhyperparasStart, FEmuVec, predictionData = NULL, fixedEffGammaAlphaBeta, errorGammaAlphaBeta,  MRAcovParasGammaAlphaBeta, maximiseOnly = FALSE, control) {
  defaultControl <- list(M = 1, randomSeed = 24,  cutForTimeSplit = 400, stepSize = 1, lowerThreshold = 3, maternCovariance = TRUE, nuggetSD = 0.00001, varyFixedEffSD = FALSE, varyMaternSmoothness = FALSE, varyErrorSD = TRUE, splitTime = FALSE, numKnotsRes0 = 20L, J = 4L, numValuesForGrid = 200, numThreads = 1, numIterOptim = 200L)
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
  if (is.null(control$lonRange)) {
    control$lonRange <- range(spacetimeData@sp@coords[ , "x"]) + c(-0.02, 0.02)
  }
  if (is.null(control$latRange)) {
    control$latRange <- range(spacetimeData@sp@coords[ , "y"]) + c(-0.02, 0.02)
  }
  if (is.null(control$timeRange)) {
    control$timeRange <- range(time(spacetimeData)) + c(-3600, 3600)
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
  # 1e10 is used as a substitute for infinity, which is not understood by the C++ code.
  if (MRAhyperparasStart$space[["smoothness"]] > 1e10) MRAhyperparasStart$space$smoothness <- 1e10
  if (MRAhyperparasStart$time[["smoothness"]] > 1e10) MRAhyperparasStart$time$smoothness <- 1e10

  for (i in names(control)) {
    defaultControl[[i]] <- control[[i]]
  }

  control <- defaultControl
  dataCoordinates <- spacetimeData@sp@coords
  predCoordinates <- predictionData@sp@coords
  predCovariates <- as.matrix(predictionData@data)

  timeRangeReshaped <- as.integer(control$timeRange)/(3600*24)
  timeBaseline <- min(timeRangeReshaped)
  timeValues <- as.integer(time(spacetimeData@time))/(3600*24) - timeBaseline # The division is to obtain values in days.
  predTime <- as.integer(time(predictionData))/(3600*24) - timeBaseline
  timeRangeReshaped <- timeRangeReshaped - timeBaseline

  covariateMatrix <- as.matrix(spacetimeData@data[, -1, drop = FALSE])
  gridPointers <- lapply(1:control$numThreads, function(x) {
    setupGridCpp(responseValues = spacetimeData@data[, 1], spCoords = dataCoordinates, predCoords = predCoordinates, obsTime = timeValues, predTime = predTime, covariateMatrix = covariateMatrix, predCovariateMatrix = predCovariates, M = control$M, lonRange = control$lonRange, latRange = control$latRange, timeRange = timeRangeReshaped, randomSeed = control$randomSeed, cutForTimeSplit = control$cutForTimeSplit, splitTime = control$splitTime, numKnotsRes0 = control$numKnotsRes0, J = control$J)$gridPointer
  })

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

  computedValues <- obtainGridValues(gridPointers = gridPointers, xStartValues = xStartValues, control = control, fixedEffSDstart = fixedEffSDstart, errorSDstart = errorSDstart, MRAhyperparasStart = MRAhyperparasStart, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, FEmuVec = FEmuVec, predictionData = predictionData, timeBaseline = timeBaseline, maximiseOnly = maximiseOnly)
  if (maximiseOnly) {
    return(computedValues$output)
  }

  hyperparaVectors <- sapply(computedValues$output, function(element) element$x)
  weightModifs <- apply(hyperparaVectors, MARGIN = 2, mvtnorm::dmvnorm, mean = computedValues$ISdistParas$mu, sigma = computedValues$ISdistParas$cov, log = TRUE)
  discreteLogJointValues <- sapply(computedValues$output, '[[', "logJointValue")
  logWeights <- discreteLogJointValues - weightModifs
  maxLogWeights <- max(logWeights)
  logPropConstantIS <- maxLogWeights + log(sum(exp(logWeights - maxLogWeights)))
  logStandardisedWeights <- logWeights - logPropConstantIS
  # Now, we obtain the marginal distribution of all mean parameters.
  print("Computing moments for marginal posterior distributions...\n")

  hyperMarginalMoments <- ComputeHyperMarginalMoments(computedValues$output, logStandardisedWeights)
  meanMarginalMoments <- ComputeMeanMarginalMoments(computedValues$output, logStandardisedWeights)
  outputList <- list(hyperMarginalMoments = hyperMarginalMoments$paraMoments, meanMarginalMoments = meanMarginalMoments, psiAndMargDistMatrix = hyperMarginalMoments$psiAndMargDistMatrix)
  print("Computing prediction moments... \n")
  if (!is.null(predictionData)) {
    outputList$predictionMoments <- ComputeKrigingMoments(computedValues$output, gridPointers[[1]], logStandardisedWeights)
    outputList$predObsOrder <- GetPredObsOrder(gridPointers[[1]]) # I picked the first one, since the grids are all copies of each other, created to ensure that there is no problem with multiple reads/writes in parallel.
  }
  print("Returning results... \n")
  outputList
}

obtainGridValues <- function(gridPointers, xStartValues, control, fixedEffSDstart, errorSDstart, MRAhyperparasStart, MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, FEmuVec, predictionData, timeBaseline, maximiseOnly = FALSE) {

  # A short optimisation first...
  storageEnvir <- new.env()
  assign(x = "x", value = NULL, envir = storageEnvir)
  assign(x = "value", value = NULL, envir = storageEnvir)
  iterCounter <- 0
  funForOptim <- function(x, envirToSaveValues) {
    iterCounter <<- iterCounter + 1
    cat("Performing iteration ", iterCounter, ".\n")
    names(x) <- names(xStartValues)
    xTrans <- exp(x)
    fixedEffArg <- fixedEffSDstart
    if (control$varyFixedEffSD) {
      fixedEffArg <- xTrans[["fixedEffSD"]]
    }
    errorArg <- errorSDstart
    if (control$varyErrorSD) {
      errorArg <- xTrans[["errorSD"]]
    }
    spSmoothnessArg <- MRAhyperparasStart$space[["smoothness"]]
    timeSmoothnessArg <- MRAhyperparasStart$time[["smoothness"]]
    if (control$varyMaternSmoothness) {
      spSmoothnessArg <- xTrans[["spSmoothness"]]
      timeSmoothnessArg <- xTrans[["timeSmoothness"]]
    }
    MRAlist <- list(space = list(rho = xTrans[["spRho"]], smoothness = spSmoothnessArg), time = list(rho = xTrans[["timeRho"]], smoothness = timeSmoothnessArg), scale = xTrans[["scale"]])
    returnedValue <- -LogJointHyperMarginal(treePointer = gridPointers[[1]], MRAhyperparas = MRAlist, fixedEffSD = fixedEffArg, errorSD = errorArg, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = control$maternCovariance, spaceNuggetSD = control$nuggetSD, timeNuggetSD = control$nuggetSD, recordFullConditional = FALSE)
    assign(x = "x", value = cbind(get(x = "x", envir = envirToSaveValues), exp(x)), envir = envirToSaveValues)
    assign(x = "value", value = c(get(x = "value", envir = envirToSaveValues), -returnedValue), envir = envirToSaveValues)
    returnedValue
  }
  gradForOptim <- function(x, envirToSaveValues) {
    names(x) <- names(xStartValues)
    numDeriv::grad(func = funForOptim, x = x, method = "simple", envirToSaveValues = envirToSaveValues)
  }
  upperBound <- rep(Inf, length(xStartValues))
  names(upperBound) <- names(xStartValues)
  upperBound <- replace(upperBound, grep(names(upperBound), pattern = "mooth"), log(50)) # This is to avoid an overflow in the computation of the Matern covariance, which for some reason does not tolerate very high smoothness values.
  upperBound <- replace(upperBound, grep(names(upperBound), pattern = "scale"), log(3000)) # This limits the scaling factor to exp(15) in the optimisation. This is to prevent computational issues in the sparse matrix inversion scheme.
  print("Optimising... \n")
  opt <- nloptr::lbfgs(x0 = log(xStartValues), lower = rep(-10, length(xStartValues)), upper = upperBound, fn = funForOptim, gr = gradForOptim, control = list(xtol_rel = 1e-3, maxeval = control$numIterOptim), envirToSaveValues = storageEnvir)
  if (!is.null(control$envirForTest)) {
    assign(x = "Hmat", value = GetHmat(gridPointers[[1]]), envir = control$envirForTest)
  }
  # cat("Optimised values:", exp(opt$par)) ;
  # stop("Stop after optimisation... \n")
  # cl <- parallel::makeForkCluster(nnodes = 4) ;
  # opt <- optimParallel::optimParallel(par = log(xStartValues), lower = rep(-10, length(xStartValues)), upper = upperBound, fn = funForOptim, gr = gradForOptim, control = list(xtol_rel = 1e-3, maxeval = control$numIterOptim), envirToSaveValues = storageEnvir)
  # opt <- nloptr::cobyla(x0 = log(xStartValues), lower = rep(-10, length(xStartValues)), upper = upperBound, fn = funForOptim, control = list(xtol_rel = 5e-4, maxeval = control$numIterOptim), envirToSaveValues = storageEnvir)
  opt$value <- -opt$value # Correcting for the inversion used to maximise instead of minimise
  sampleWeights <- exp(storageEnvir$value - max(storageEnvir$value))
  sampleWeights <- sampleWeights/sum(sampleWeights)
  varCovar <- cov.wt(x = t(storageEnvir$x), wt = sampleWeights)$cov
  solution <- exp(opt$par)
  if (maximiseOnly) {
    names(solution) <- names(xStartValues)
    return(solution)
  }

  if (length(gridPointers) > 1) {
    require(doParallel)
    no_cores <- length(gridPointers)
    doParallel::registerDoParallel(cores = no_cores)# Shows the number of Parallel Workers to be used
  }
  gridFct <- function() {
    smallMVR <- function() {
      container <- NULL
      repeat {
        container <- drop(mvtnorm::rmvnorm(n = 1, mean = exp(opt$par), sigma = varCovar))
        if (all(container > 0)) break
      }
      container
    }
    paraGrid <- t(replicate(n = control$numValuesForGrid, expr = smallMVR()))
    colnames(paraGrid) <- names(xStartValues)
    groupAssignments <- rep(1:length(gridPointers), length.out = nrow(paraGrid))
    listForParallel <- lapply(1:min(length(gridPointers), nrow(paraGrid)), function(clIndex) {
      indicesForNode <- which(groupAssignments == clIndex)
      list(paraGrid = paraGrid[indicesForNode,, drop = FALSE])
    })
    output <- vector("list", length = nrow(paraGrid))
    if (length(gridPointers) == 1) {
      # output <- lapply(1:nrow(paraGrid), funForGridEst, paraGrid = paraGrid, treePointer = gridPointers[[1]], MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, fixedEffSDstart = fixedEffSDstart, errorSDstart = errorSDstart, MRAhyperparasStart = MRAhyperparasStart, FEmuVec = FEmuVec, predictionData = predictionData, timeBaseline = timeBaseline, computePrediction = TRUE, control = control)
      for (i in 1:length(output)) {
        cat("Processing point ", i, ".\n")
        output[[i]] <- funForGridEst(index = i, paraGrid = paraGrid, treePointer = gridPointers[[1]], MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, fixedEffSDstart = fixedEffSDstart, errorSDstart = errorSDstart, MRAhyperparasStart = MRAhyperparasStart, FEmuVec = FEmuVec, predictionData = predictionData, timeBaseline = timeBaseline, computePrediction = TRUE, control = control)
      }
    } else {
      output <- foreach::foreach(var1 = seq_along(listForParallel), .inorder = FALSE) %dopar% {
        lapply(1:nrow(listForParallel[[var1]]$paraGrid), funForGridEst, paraGrid = listForParallel[[var1]]$paraGrid, treePointer = gridPointers[[var1]], MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, fixedEffSDstart = fixedEffSDstart, errorSDstart = errorSDstart, MRAhyperparasStart = MRAhyperparasStart, FEmuVec = FEmuVec, predictionData = predictionData, timeBaseline = timeBaseline, computePrediction = TRUE, control = control)
      }
      output <- do.call("c", output)
    }
    # list(output = output, optimPoints = list(x = storageEnvir$x, value = storageEnvir$value))
    list(output = output)
  }
  print("Computing values on the grid...")
  valuesOnGrid <- gridFct()
  print("Grid complete... \n")
  if (length(gridPointers) > 1) {
    doParallel::stopImplicitCluster()
  }
  keepIndices <- sapply(valuesOnGrid$output, function(x) class(x$logJointValue) == "numeric")
  valuesOnGrid$output <- valuesOnGrid$output[keepIndices]
  valuesOnGrid$ISdistParas <- list(mu = exp(opt$par), cov = varCovar)
  valuesOnGrid
}

funForGridEst <- function(index, paraGrid, treePointer, predictionData, MRAcovParasGammaAlphaBeta, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, fixedEffSDstart, errorSDstart, MRAhyperparasStart, FEmuVec, timeBaseline, computePrediction, control) {
  if (is.null(control$batchSizePredict)) {
    control$batchSizePredict <- 500
  }
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
  print("Computing log-joint... \n")
  logJointValue <- tryCatch(expr = LogJointHyperMarginal(treePointer = treePointer, MRAhyperparas = MRAlist, fixedEffSD = fixedEffArg, errorSD = errorArg, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta = fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = control$maternCovariance, spaceNuggetSD = control$nuggetSD, timeNuggetSD = control$nuggetSD, recordFullConditional = FALSE), error = function(e) e)
  print("Done computing log-joint! \n")
  errorSD <- ifelse(any(grepl(pattern = "error", x = names(x))), x[[grep(pattern = "error", x = names(x), value = TRUE)]], errorSDstart)
  fixedEffSD <- ifelse(any(grepl(pattern = "fixedEff", x = names(x))), x[[grep(pattern = "fixedEff", x = names(x), value = TRUE)]], fixedEffSDstart)
  aList <- list(x = x, errorSD = errorSD, fixedEffSD = fixedEffSD, MRAhyperparas = MRAlist, logJointValue = logJointValue)
  if (!is.null(predictionData) & computePrediction) {
    timeValues <- as.integer(time(predictionData@time))/(3600*24) - timeBaseline # The division is to obtain values in days.
    print("Computing statistics for predictions... \n")
    aList$CondPredStats <- ComputeCondPredStats(treePointer, predictionData@sp@coords, timeValues, as.matrix(predictionData@data), batchSize = control$batchSizePredict)
    print("Done computing statistics for predictions! \n")
  }
  # Running LogJointHyperMarginal stores in the tree pointed by gridPointer the full conditional mean and SDs when recordFullConditional = TRUE. We can get them with the simple functions I call now.
  aList$FullCondMean <- GetFullCondMean(treePointer)
  aList$FullCondSDs <- GetFullCondSDs(treePointer)
  aList
}

ComputeHyperMarginalMoments <- function(hyperparaList, logISmodProbWeights) {
  domainCheck <- sapply(hyperparaList, function(x) x$logJointValue > -Inf)
  hyperparaList <- hyperparaList[domainCheck]
  psiAndMargDistMatrix <- t(sapply(seq_along(hyperparaList), function(hyperparaIndex) c(unlist(hyperparaList[[hyperparaIndex]]$MRAhyperparas), fixedEffSD = hyperparaList[[hyperparaIndex]]$fixedEffSD, errorSD = hyperparaList[[hyperparaIndex]]$errorSD, jointValue = exp(logISmodProbWeights[hyperparaIndex]))))
  rownames(psiAndMargDistMatrix) <- NULL
  funToGetParaMoments <- function(hyperparaIndex) {
    meanValue <- sum(psiAndMargDistMatrix[, hyperparaIndex] * psiAndMargDistMatrix[, ncol(psiAndMargDistMatrix)])
    sdValue <- sqrt(sum(psiAndMargDistMatrix[, hyperparaIndex]^2 * psiAndMargDistMatrix[, ncol(psiAndMargDistMatrix)]) - meanValue^2)
    c(mean = meanValue, StdDev = sdValue)
  }
  paraMoments <- t(sapply(1:(ncol(psiAndMargDistMatrix) - 1), FUN = funToGetParaMoments))

  colnames(paraMoments) <- c("Mean", "StdDev")
  rownames(paraMoments) <- head(colnames(psiAndMargDistMatrix), n = -1)
  list(paraMoments = as.data.frame(paraMoments), psiAndMargDistMatrix = psiAndMargDistMatrix)
}

ComputeMeanMarginalMoments <- function(hyperparaList, logISmodProbWeights) {

  numMeanParas <- length(hyperparaList[[1]]$FullCondMean)

  marginalMeans <- sapply(1:numMeanParas, function(paraIndex) {
    meanVector <- sapply(hyperparaList, function(x) x$FullCondMean[[paraIndex]])
    sum(meanVector * exp(logISmodProbWeights))
  })
  marginalSecondMoments <- sapply(1:numMeanParas, function(paraIndex) {
    meanVector <- sapply(hyperparaList, function(x) x$FullCondMean[[paraIndex]])
    sdVector <- sapply(hyperparaList, function(x) x$FullCondSDs[[paraIndex]])
    secondMomentVec <- sdVector^2 + meanVector^2
    sum(secondMomentVec * exp(logISmodProbWeights))
  })
  marginalSDs <- marginalSecondMoments - marginalMeans^2
  data.frame(Mean = marginalMeans, StdDev = marginalSDs)
}

ComputeKrigingMoments <- function(hyperparaList, treePointer, logISmodProbWeights) {

  termsForMean <- lapply(seq_along(hyperparaList), function(index) {
    drop(hyperparaList[[index]]$CondPredStats$Hmean * exp(logISmodProbWeights[[index]]))
  })
  krigingMeans <- Reduce("+", termsForMean)
  termsForVarE1 <- lapply(seq_along(hyperparaList), FUN = function(index) {
    drop(hyperparaList[[index]]$CondPredStats$Hmean^2 * exp(logISmodProbWeights[[index]]))
  })
  varE <- Reduce('+', termsForVarE1) - krigingMeans^2

  termsForEvar <- lapply(seq_along(hyperparaList), function(index) {
    drop(hyperparaList[[index]]$CondPredStats$Evar * exp(logISmodProbWeights[[index]]))
  })
  Evar <- Reduce("+", termsForEvar)
  predObsOrder <- GetPredObsOrder(treePointer = treePointer)
  list(predictMeans = krigingMeans[order(predObsOrder)], predictSDs = sqrt(varE + Evar)[order(predObsOrder)])
}

covFunctionBiMatern <- function(rangeParaSpace = 10, rangeParaTime = 10) {
  function(spacetime1, spacetime2) {
    euclidDist <- sp::spDists(spacetime1@sp, spacetime2@sp)
    timeDist <- outer(zoo::index(spacetime2@time), zoo::index(spacetime1@time), function(x, y) as.numeric(abs(difftime(x, y, units = "days"))))
    fields::Exponential(euclidDist, range = 10)*t(fields::Exponential(timeDist, range = 10))
  }
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

  logDetFun <- function(cholM) {
    cholM <- as(cholM, "dgCMatrix")
    2*sum(log(diag(cholM)))
  }

  LogJointHyperMarginalToWrap(treePointer = treePointer, MRAhyperparas = MRAhyperparas, fixedEffSD = fixedEffSD, errorSD = errorSD, MRAcovParasGammaAlphaBeta = MRAcovParasGammaAlphaBeta, FEmuVec = FEmuVec, fixedEffGammaAlphaBeta =  fixedEffGammaAlphaBeta, errorGammaAlphaBeta = errorGammaAlphaBeta, matern = matern, spaceNuggetSD = spaceNuggetSD, timeNuggetSD = timeNuggetSD, recordFullConditional = TRUE)
}
