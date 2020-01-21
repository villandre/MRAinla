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

MRA_INLA <- function(spacetimeData, hyperStart, fixedHyperValues, hyperpriorPars, FEmuVec, predictionData = NULL,  maximiseOnly = FALSE, control) {
  # CHECKS #########################
  if (!is.null(predictionData)) {
    if (!identical(colnames(spacetimeData@data)[-1], colnames(predictionData@data))) {
      stop("Mismatch between covariates in training and test data. \n")
    }
  }
  ##################################
  # DEFINING CONTROL PARA.##########
  defaultControl <- list(Mlon = 1, Mlat = 1, Mtime = 1, randomSeed = 24,  cutForTimeSplit = 400, nuggetSD = 0.00001, splitTime = FALSE, numKnotsRes0 = 20L, J = 4L, numValuesForIS = 200, numIterOptim = 200L, distMethod = "haversine", normalHyperprior = FALSE)
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

  gridPointer <- setupGridCpp(responseValues = spacetimeData@data[, 1], spCoords = dataCoordinates, predCoords = predCoordinates, obsTime = timeValues, predTime = predTime, covariateMatrix = covariateMatrix, predCovariateMatrix = predCovariates, Mlon = control$Mlon, Mlat = control$Mlat, Mtime = control$Mtime, lonRange = control$lonRange, latRange = control$latRange, timeRange = timeRangeReshaped, randomSeed = control$randomSeed, cutForTimeSplit = control$cutForTimeSplit, splitTime = control$splitTime, numKnotsRes0 = control$numKnotsRes0, J = control$J, distMethod = control$distMethod, MaternParsHyperpars = hyperpriorPars[c("space", "time", "scale")], fixedEffParsHyperpars = hyperpriorPars$fixedEffSD, errorParsHyperpars = hyperpriorPars$errorSD, FEmuVec = FEmuVec, nuggetSD = control$nuggetSD, normalHyperprior = control$normalHyperprior)$gridPointer

  # First we compute values relating to the hyperprior marginal distribution...

  computedValues <- .obtainGridValues(gridPointer = gridPointer, hyperStart = hyperStart, fixedHyperValues = fixedHyperValues, predictionData = predictionData, timeBaseline = timeBaseline, maximiseOnly = maximiseOnly, control = control)
  if (maximiseOnly) {
    return(computedValues$output)
  }

  hyperparaVectors <- sapply(computedValues$output, function(element) element$x)
  weightModifs <- apply(hyperparaVectors, MARGIN = 2, mvtnorm::dmvnorm, mean = computedValues$ISdistParas$mu, sigma = computedValues$ISdistParas$cov, log = TRUE)
  discreteLogJointValues <- sapply(computedValues$output, '[[', "logJointValue")
  logWeights <- discreteLogJointValues - weightModifs - log(length(discreteLogJointValues))
  maxLogWeights <- max(logWeights)
  logPropConstantIS <- maxLogWeights + log(sum(exp(logWeights - maxLogWeights)))
  logStandardisedWeights <- logWeights - logPropConstantIS
  # Now, we obtain the marginal distribution of all mean parameters.
  cat("Computing moments for marginal posterior distributions...\n")

  hyperMarginalMoments <- ComputeHyperMarginalMoments(computedValues$output, logStandardisedWeights)
  meanMarginalMoments <- ComputeMeanMarginalMoments(computedValues$output, logStandardisedWeights)
  outputList <- list(hyperMarginalMoments = hyperMarginalMoments$paraMoments, meanMarginalMoments = meanMarginalMoments, psiAndMargDistMatrix = hyperMarginalMoments$psiAndMargDistMatrix)
  cat("Computing prediction moments... \n")
  if (!is.null(predictionData)) {
    outputList$predictionMoments <- ComputeKrigingMoments(computedValues$output, gridPointer, logStandardisedWeights, control = control)
    outputList$predObsOrder <- GetPredObsOrder(gridPointer) # I picked the first one, since the grids are all copies of each other, created to ensure that there is no problem with multiple reads/writes in parallel.
    if (identical(control$IScompleted, TRUE)) { # m_predObsOrder needs to be obtained from the hard drive, as CreateHmatrixPred was not run.
      outputList$predObsOrder <- get(load(paste(control$folderToSaveISpoints, "/predObsOrder.Rdata", sep = "")))
    }
  }
  cat("Returning results... \n")
  outputList
}

.obtainGridValues <- function(gridPointer, hyperStart, fixedHyperValues, predictionData, timeBaseline, maximiseOnly = FALSE, control) {
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
    returnedValue <- -LogJointHyperMarginal(treePointer = gridPointer, hyperparaValues = hyperList, recordFullConditional = FALSE, processPredictions = FALSE)
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
    varCovar <- solve(precisionMat)
    if (!is.null(control$fileToSaveOptOutput)) {
      save(opt, file = control$fileToSaveOptOutput)
      filenameForVarCovar <- paste(substr(control$fileToSaveOptOutput, start = 1, stop = gregexpr(pattern = ".Rdata", text = control$fileToSaveOptOutput)[[1]] - 1), "_ISvarCovar.Rdata", sep = "")
      save(varCovar, file = filenameForVarCovar)
    }
    cat("Optimised values:", solution)
  } else {
    load(control$fileToSaveOptOutput)
    solution <- opt$par
    if (!control$normalHyperprior) {
      solution <- exp(opt$par)
    }
    varCovarFilename <- paste(substr(control$fileToSaveOptOutput, start = 1, stop = gregexpr(pattern = ".Rdata", text = control$fileToSaveOptOutput)[[1]] - 1), "_ISvarCovar.Rdata", sep = "")
    load(varCovarFilename) # Restores varCovar
  }

  if (!is.null(control$envirForTest)) {
    assign(x = "Hmat", value = GetHmat(gridPointer), envir = control$envirForTest)
  }

  opt$value <- -opt$value # Correcting for the inversion used to maximise instead of minimise

  if (maximiseOnly) {
    names(solution) <- names(xStartValues)
    return(solution)
  }

  cat("Running IS algorithm... \n")
  valuesOnGrid <- .gridFct(distMode = solution, ISvarCovar = varCovar, gridPointer = gridPointer, namesXstartValues = names(xStartValues), fixedHyperValues = fixedHyperValues, control = control)
  cat("IS algorithm completed... \n")

  keepIndices <- sapply(valuesOnGrid$output, function(x) class(x$logJointValue) == "numeric")
  valuesOnGrid$output <- valuesOnGrid$output[keepIndices]
  valuesOnGrid$ISdistParas <- list(mu = solution, cov = varCovar)
  valuesOnGrid
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

.gridFct <- function(distMode, ISvarCovar, gridPointer, namesXstartValues, fixedHyperValues, control) {
  if (!control$normalHyperprior) {
    smallMVR <- function() {
      container <- NULL
      repeat {
        container <- drop(mvtnorm::rmvnorm(n = 1, mean = distMode, sigma = ISvarCovar))
        if (all(container > 0)) break
      }
      container
    }
    paraGrid <- t(replicate(n = control$numValuesForIS, expr = smallMVR()))
  } else {
    paraGrid <- mvtnorm::rmvnorm(n = control$numValuesForIS, mean = distMode, sigma = ISvarCovar)
  }
  colnames(paraGrid) <- namesXstartValues
  output <- vector("list", length = nrow(paraGrid))
  startAtIter <- 1
  if (tryCatch(dir.exists(control$folderToSaveISpoints), error = function(e) FALSE)) {
    cat("Loading previously processed IS points... \n")
    loadResult <- lapply(list.files(control$folderToSaveISpoints, full.names = TRUE, pattern = "ISoutput"), function(x) get(load(x)))
    output[1:length(loadResult)] <- loadResult
    if (identical(control$IScompleted, TRUE)) { # If we indicate that the IS is completed, the function will simply stop sampling points and use saved points only.
      output <- output[1:length(loadResult)]
    }
    startAtIter <- length(loadResult) + 1
  }
  if (startAtIter <= length(output)) {
    for (i in startAtIter:length(output)) {
      cat("Processing grid value ", i, "... \n")
      xVec <- paraGrid[i, ]
      names(xVec) <- colnames(paraGrid)
      if (control$normalHyperprior) {
        xVec <- exp(xVec)
      }
      output[[i]] <- .funForGridEst(xNonLogScale = xVec, treePointer = gridPointer, fixedHyperValues = fixedHyperValues, computePrediction = TRUE, control = control)
      if (!is.null(control$folderToSaveISpoints)) {
        if (!dir.exists(control$folderToSaveISpoints)) {
          dir.create(path = control$folderToSaveISpoints)
        }
        if (i == 1) { # On the first iteration, the vector to restore the order of predictions must be saved to allow for the predicted values to be re-ordered after a resume in which CreateHmatrixPred is not called. This situation occurs when we specify control$IScompleted=TRUE.
          predOrder <- GetPredObsOrder(gridPointer)
          filenameForPredOrder <- paste(control$folderToSaveISpoint, "/predObsOrder.Rdata", sep = "")
          save(predOrder, file = filenameForPredOrder, compress = TRUE)
        }
        objectToSave <- output[[i]]
        filename <- paste(control$folderToSaveISpoints, "/", "ISoutput", i, ".Rdata", sep = "")
        save(objectToSave, file = filename, compress = TRUE)
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

  logJointValue <- tryCatch(expr = LogJointHyperMarginal(treePointer = treePointer, hyperparaValues = hyperList, recordFullConditional = FALSE, processPredictions = TRUE), error = function(e) e)
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
  # Running LogJointHyperMarginal stores in the tree pointed by gridPointer the full conditional mean and SDs when recordFullConditional = TRUE. We can get them with the simple functions I call now.
  aList$FullCondMean <- GetFullCondMean(treePointer)
  aList$FullCondSDs <- GetFullCondSDs(treePointer)
  aList
}

ComputeHyperMarginalMoments <- function(hyperparaList, logISmodProbWeights) {
  domainCheck <- sapply(hyperparaList, function(x) x$logJointValue > -Inf)
  hyperparaList <- hyperparaList[domainCheck]
  psiAndMargDistMatrix <- t(sapply(seq_along(hyperparaList), function(hyperparaIndex) c(unlist(hyperparaList[[hyperparaIndex]]$MaternHyperpars), fixedEffSD = hyperparaList[[hyperparaIndex]]$fixedEffSD, errorSD = hyperparaList[[hyperparaIndex]]$errorSD, jointValue = exp(logISmodProbWeights[hyperparaIndex]))))
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
  marginalSDs <- sqrt(marginalSecondMoments - marginalMeans^2)
  data.frame(Mean = marginalMeans, StdDev = marginalSDs)
}

ComputeKrigingMoments <- function(hyperparaList, treePointer, logISmodProbWeights, control) {

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
  if (identical(control$IScompleted, TRUE)) {
    predObsOrder <- get(load(paste(control$folderToSaveISpoints, "/predObsOrder.Rdata", sep = "")))
  }
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

LogJointHyperMarginal <- function(treePointer, hyperparaValues, recordFullConditional, processPredictions = FALSE) {
  LogJointHyperMarginalToWrap(treePointer = treePointer, MaternHyperpars = hyperparaValues[c("space", "time", "scale")], fixedEffSD = hyperparaValues$fixedEffSD, errorSD = hyperparaValues$errorSD, recordFullConditional = TRUE, processPredictions = processPredictions)
}
