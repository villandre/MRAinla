# spaceDistFun must take as an argument a matrix or data.frame with a column named latitude and another called longitude.

.SimulateSpacetimeData <- function(numObsPerTimeSlice = 225, covFunction, lonRange, latRange,  timeValuesInPOSIXct, covariateGenerationFctList, errorSD, spaceDistFun, timeDistFun, FEvalues) {
  numSlotsPerRow <- ceiling(sqrt(numObsPerTimeSlice))
  slotCoordinates <- sapply(list(longitude = lonRange, latitude = latRange), FUN = function(x) {
    width <- abs(diff(x)/numSlotsPerRow)
    seq(from = min(x) + width/2, to = max(x) - width/2, length.out = numSlotsPerRow)
  }, USE.NAMES = TRUE)

  allSpaceCoordinates <- as.data.frame(expand.grid(as.data.frame(slotCoordinates)))
  numToRemove <- nrow(allSpaceCoordinates) - numObsPerTimeSlice

  if (numToRemove > 0) {
    obsToRemove <- (nrow(allSpaceCoordinates) - numToRemove + 1):nrow(allSpaceCoordinates)
    allSpaceCoordinates <- allSpaceCoordinates[-obsToRemove, ]
  }

  coordinates <- allSpaceCoordinates[rep(1:nrow(allSpaceCoordinates), length(timeValuesInPOSIXct)), ]
  coordinates$time <- rep(timeValuesInPOSIXct, each = numObsPerTimeSlice)

  covariateMatrix <- cbind(1, do.call("cbind", lapply(covariateGenerationFctList, function(x) {
    covariateValues <- x(coordinates[, c("longitude", "latitude")], coordinates[, "time"])
    if (is.null(dim(covariateValues))) {
      dim(covariateValues) <- c(length(covariateValues), 1)
    }
    covariateValues
  })))

  spatialDistMatrix <- spaceDistFun(coordinates[ , c("longitude", "latitude")])
  timeDistMatrix <- timeDistFun(coordinates[, "time"])
  covarianceMat <- covFunction(spatialDistMatrix, timeDistMatrix)
  meanVector <- drop(covariateMatrix %*% FEvalues)
  fieldValues <- drop(mvtnorm::rmvnorm(n = 1, mean = meanVector, sigma = covarianceMat)) + rnorm(n = length(meanVector), mean = 0, sd = errorSD)
  dataForObject <- cbind(y = fieldValues, as.data.frame(covariateMatrix[, -1]))
  colnames(dataForObject) <- c("y", paste("Covariate", 1:(length(FEvalues) - 1), sep = ""))
  spacetimeObj <- spacetime::STIDF(sp = sp::SpatialPoints(coordinates[, c("longitude", "latitude")]), data = dataForObject, time = coordinates$time)
  colnames(spacetimeObj@sp@coords) <- c("longitude", "latitude")
  spacetimeObj
}

.plotSpacetimeData <- function(spacetimeData, fontsize) {

  width <- ceiling(sqrt(nrow(spacetimeData@sp@coords)/length(unique(time(spacetimeData@time)))))
  padding <- 0.01
  landRaster <- raster::raster(nrows = width, ncols = width, xmn = min(spacetimeData@sp@coords[ , 1]) - padding, xmx = max(spacetimeData@sp@coords[ , 1]) + padding, ymn = min(spacetimeData@sp@coords[ , 2]) - padding, ymx = max(spacetimeData@sp@coords[ , 2]) + padding)
  getTestRaster <- function(dayIndex) {

    desiredTime <- sort(unique(time(spacetimeData@time)))[[dayIndex]]

    timeIndicesTest <- which(time(spacetimeData@time) == desiredTime)

    subSpacetimeData <- spacetimeData[timeIndicesTest]

    testPoints <- sp::SpatialPoints(subSpacetimeData@sp@coords)

    raster::rasterize(x = testPoints, y = landRaster, field = subSpacetimeData@data$y)
  }
  rasterList <- lapply(seq_along(unique(time(spacetimeData@time))), getTestRaster)
  stackedRasters <- raster::stack(rasterList)
  raster::spplot(stackedRasters, scales = list(draw = FALSE),
         xlab = "Longitude", ylab = "Latitude",
         names.attr = as.character(unique(time(spacetimeData))),
         par.settings = list(fontsize=list(text=fontsize)))
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
