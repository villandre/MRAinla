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

MRA_INLA <- function(spaceTimeList, spaceTimeCovFct, M, gridRasterList, numKnots, hyperpriorFunList) {

}

.logLikFun <- function(spaceTimeList, meanValueList, sigmaErrorValue) {
  sum(mapply(grid = spaceTimeList, meanValues = meanValueList, FUN = function(grid, meanValues) {
    sum(dnorm(x = grid@data, mean = meanValues, sd = sigmaErrorValue, log = TRUE))
  }))
}

.logGMRFprior <- function(paraValuesList, SigmaValuesList) {
  sum(mapply(paraValues = paraValuesList, varPar = SigmaValuesList, FUN = function(paraValues, varPar) {
    if (!is.matrix(varPar)) {
      varPar <- varPar*diag(length(paraValues))
    }
    mvtnorm::dmvnorm(x = paraValues, sigma = varPar, log = TRUE)
  }))
}

.createPriorFunction <- function(gridList, knotsList, covFct, resolution = length(gridList)) {
  function(spaceTime1, spaceTime2) {
    currentVfun <- covFct
    currentKfun <- function() {
      covFct(knotsList[[1]], knotsList[[1]])
    }
    currentBfun <- function(spaceTime) {
      covFct(spaceTime, knotsList[[1]])
    }
    initValues <- list(v = currentVfun(spaceTime1, spaceTime2), K = currentKfun(), bList = list(currentBfun(spaceTime1), currentBfun(spaceTime2)))
    if (identical(resolution, 0)) {
      return(initValues)
    }
    funForApply <- function(resIndex) {
      if (!.sameGridSection(spaceTime1, spaceTime2, gridList[[resIndex]])) {
        return(list(v = 0, K = 0, bList = 0))
      }
      knotsInRegion <- .getPointsInRegion(gridList[[resIndex]], spaceTime1, knotsList[[resIndex+1]])
      newVfun <- .setupVrecursionStep(spacetimegridObj = gridList[[resIndex]], vFun = currentVfun, bFun = currentBfun, Kfun = currentKfun)
      newBfun <- function(spaceTime1) {
        .getBvec(spaceTimeCoord = spaceTime1, knotPositions = knotsInRegion, vFun = newVfun)
      }
      newKfun <- function() {
        .getKmat(knotPositions = knotsInRegion, vFun = newVfun)
      }
      returnValues <- list(v = newVfun(spaceTime1,spaceTime2), K = newKfun(), bList = list(newBfun(spaceTime1), newBfun(spaceTime2)))
      currentVfun <<- newVfun
      currentBfun <<- newBfun
      currentKfun <<- newKfun
      return(returnValues)
    }
    fittedValues <- lapply(1:resolution, FUN = funForApply)
    c(list(initValues), fittedValues)
  }
}

.setupVrecursionStep <- function(spacetimegridObj, vFun, bFun, Kfun) {
  function(spaceTime1, spaceTime2) {
    if(!.sameGridSection(spaceTime1, spaceTime2, spacetimegridObj)) {
      return(0)
    }
    vFun(spaceTime1, spaceTime2) - t(bFun(spaceTime1))%*%Kfun()%*%bFun(spaceTime2)
  }
}

.initVfun <- function(knots, covFct) {
  function(spaceTime1, spaceTime2) {
    bFun <- function(spaceTime) {
      covFct(spaceTime, knots)
    }
    Kfun <- function()
    {
      Matrix::chol2inv(Matrix::chol(covFct(knots, knots)))
    }
    list(v = covFct(spaceTime1, spaceTime2), K = Kfun(), bValues = list(bFun(spaceTime1), bFun(spaceTime2)))
  }
}

.getBvec <- function(spaceTimeCoord, knotPositions, vFun) {
  sapply(1:Npoints(knotPositions), FUN = function(knotIndex) {
    vFun(spaceTimeCoord, knotPositions[knotIndex])
  })
}

.getKmat <- function(knotPositions, vFun) {
  KmatrixFinal <- matrix(0, nrow(knotPositions), nrow(knotPositions))
  # Handling the off diagonal elements
  mapply(rowIndex = row(diag(nrow(knotPositions)))[lower.tri(diag(nrow(knotPositions)))], colIndex = col(diag(nrow(knotPositions)))[lower.tri(diag(nrow(knotPositions)))], FUN =function(rowIndex, colIndex) {
    KmatrixFinal[rowIndex, colIndex] <<- vFun(knotPositions[rowIndex], knotPositions[colIndex])
    invisible(NULL)
  }, SIMPLIFY = FALSE)
  KmatrixFinal <- KmatrixFinal + t(KmatrixFinal)
  diag(KmatrixFinal) <- sapply(1:nrow(knotPositions), FUN = function(x) {
    vFun(knotPositions[x], knotPositions[x])
  })
  KmatrixFinal
}

#' Constructs an irregular grid on a specified spatiotemporal range.
#'
#' Constructing a spatiotemporal grid will be necessary prior to fitting the spatiotemporal MRA
#' proposed by Villandre et al.
#'
#' @param lonBreaks vector indicating where the longitude boudaries are located
#' @param latBreaks vector indicating where the latitude boudaries are located
#' @param timeBreaks vector indicating where the temporal boudaries are located
#' @param numKnots either the number of knots in all elements of the grids, or a vector of size equal to
#' the length of gridRasterList giving the number of knots for each resolution, or a list of vectors giving
#' the number of knots in each region of the grid for each resolution
#' @param hyperPriorFunList named list of functions with one argument specifying the hyperprior distributions
#'
#' @details A 3D irregular grid can be represented parsemoniously by a list of 3 vectors, one #' for each dimension. Under this scheme, spatiotemporal regions at a given resolution are
#' like bricks. This is a constructor function.
#'
#' @return A Spacetimegrid object.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

gridConstructor <- function(lonBreaks, latBreaks, timeBreaks, lonExtent, latExtent, timeExtent) {
  grids <- mapply(breaks = list(lonBreaks, latBreaks, timeBreaks), extent = list(lonExtent, latExtent, timeExtent), FUN = function(breaks, extent) {
    list(breaks = breaks, extent = extent)
  }, SIMPLIFY = FALSE)
  names(grids) <- c("longitude", "latitude", "time")
  class(grids) <- "Spacetimegrid"
  grids
}

#' A custom print function for the Spacetimegrid class.
#'
#' Prints a Spacetimegrid object.
#'
#' @param x Spacetimegrid object
#'
#' @details
#'
#' @return Print function, doesn't return anything.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

print.Spacetimegrid <- function(x) {
  cat("Spacetime grid object \n")
  cat("Longitude extent: ", x$longitude$extent , "\n")
  cat("Latitude extent: ", x$latitude$extent, "\n")
  cat("Time extent: ")
  cat(paste(x$time$extent, collapse = " "), "\n \n")
  cat("Grid dimensions: (", length(x$longitude$breaks) + 1, " ", length(x$latitude$breaks)+1, " ", length(x$time$breaks) + 1, ")\n")
  cat("Note: This gives the number of regions, which may have different volumes.\n\n")
  invisible(NULL)
}

#' A custom plot function for the Spacetimegrid class.
#'
#' Plots the spatial grid within the Spacetimegrid object, plus, optionally, observations and knots.
#'
#' @param x Spacetimegrid object
#' @param observationsAsPoints SpatialPoints or SpatialPointsDataFrame object containing the
#' observations
#' @param knotsAsPoints SpatialPoints or SpatialPointsDataFrame object containing the
#' knots
#'
#' @details The function plots a 2D projection of the grid, flattening the dataset with respect to time.
#'
#' @return Plot function, doesn't return anything worthwhile.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export

plot.Spacetimegrid <- function(x, observationsAsPoints, knotsAsPoints) {
  combinedSP <- NULL
  plotColours <- NULL
  if (!is.null(observationsAsPoints)) {
    observSP <- SpatialPointsDataFrame(observationsAsPoints, data = data.frame(identity = rep("red", nrow(observationsAsPoints@coords))))
    combinedSP <- observSP
    plotColours <- combinedSP@data$identity
  }
  if (!is.null(knotsAsPoints)) {
    knotSP <- SpatialPointsDataFrame(knotsAsPoints, data = data.frame(identity = rep("green", nrow(knotsAsPoints@coords))))
    combinedSP <- raster::bind(combinedSP, knotSP)
    plotColours <- combinedSP@data$identity
  }
  plot(combinedSP, col = as.character(plotColours), xlab = "Longitude", ylab = "Latitude", xlim = x$longitude$extent, ylim = x$latitude$extent, cex = 0.8, pch = 18)
  abline(v = c(x$longitude$extent[1], x$longitude$breaks, x$longitude$extent[2]), lty = 3, col = "blue4")
  abline(h = c(x$latitude$extent[1], x$latitude$breaks, x$latitude$extent[2]), lty = 3, col = "blue4")
}

# observations are coded as SpatialPointsDataFrame with 3D coordinates, longitude-latitude-time
# In each dimension, we have ranges defined [.,.), i.e. closed on the left, open on the right.

.getPointsInRegion <- function(sptGrid, spacetimePoint, observations) {
  lapply(c("longitude", "latitude", "time"), FUN = function(dimensionName) {
    upperPos <- match(TRUE, sptGrid[[dimensionName]]$breaks > .getSpacetimeDim(spacetimePoint, dimension = dimensionName))
    coorRange <- c(sptGrid[[dimensionName]]$breaks[upperPos - 1], sptGrid[[dimensionName]]$breaks[upperPos])
    observations <<- observations[(.getSpacetimeDim(observations, dimension =  dimensionName) >= coorRange[[1]]) & (.getSpacetimeDim(observations, dimension =  dimensionName) < coorRange[[2]])]
    invisible(NULL)
  })
  observations
}

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

#' Add breaks in space-time grid.
#'
#' This function is used to refine a resolution grid.
#'
#' @param spacetimegridObj Spacetimegrid object
#' @param dimension the name of the dimension where breaks will be added, either "longitude",
#' "latitude", or "time"
#' @param breaks vector indicating where the new breakpoints should be located.
#'
#' @details This function is meant to create an embedded resolution.
#'
#' @return A Spacetimegrid object with the new breaks.
#'
#' @examples
#' \dontrun{
#' INPUT_AN_EXAMPLE()
#' }
#' @export
addBreaks <- function(spacetimegridObj, dimension = c("longitude", "latitude", "time"), breaks) {
  dimension <- dimension[[1]]
  spacetimegridObj[[dimension]]$breaks <- sort(c(spacetimegridObj[[dimension]]$breaks, breaks))
  spacetimegridObj
}

.sameGridSection <- function(x, y, spacetimegridObj) {
  sameSection <- FALSE
  testResultTime <- identical(match(TRUE, spacetimegridObj[["time"]]$breaks > index(x)), match(TRUE, spacetimegridObj[["time"]]$breaks > index(y)))
  if (testResultTime) {
    testResultsSpace <- sapply(c("longitude", "latitude"), FUN = function(dimensionName) {
      xUpperPos <- match(TRUE, spacetimegridObj[[dimensionName]]$breaks > .getSpacetimeDim(x, dimensionName))
      yUpperPos <- match(TRUE, spacetimegridObj[[dimensionName]]$breaks > .getSpacetimeDim(y, dimensionName))
      identical(xUpperPos, yUpperPos)
    })
    if (all(testResultsSpace)) {
      sameSection <- TRUE
    }
  }
  sameSection
}

spacetimeListConvertToPoints <- function(valuesList, timeValues=NULL, regular = FALSE) {
  if (regular) {
    if (class(valuesList[[1]]) == "SpatialPoints") {
      return(spacetime::STF(sp = valuesList[[1]], time = timeValues))
    }
    valuesFrames <- lapply(valuesList, FUN = function(x) x@data)
    dataReformat <- do.call("rbind", valuesFrames)
    return(spacetime::STFDF(valuesList[[1]], time = timeValues, data = dataReformat))
  }
  allPoints <- do.call(what = raster::bind, valuesList)
  lengthVec <- sapply(valuesList, length)
  extendedTimeValues <- rep(timeValues, lengthVec)

  if (class(allPoints) == "SpatialPoints") {
    return(spacetime::STI(sp = allPoints, time = extendedTimeValues))
  }
  spacetime::STIDF(sp = allPoints, time = extendedTimeValues, data = allPoint@data)
}

Npoints <- function(spacetimeObj) {
  nrow(spacetimeObj@sp@coords)
}
