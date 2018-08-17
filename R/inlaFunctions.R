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

logLikFun <- function(spaceTimeList, meanValueList, sigmaErrorValue) {
  sum(mapply(grid = spaceTimeList, meanValues = meanValueList, FUN = function(grid, meanValues) {
    sum(dnorm(x = grid@data, mean = meanValues, sd = sigmaErrorValue, log = TRUE))
  }))
}

logGMRFprior <- function(paraValuesList, SigmaValuesList) {
  sum(mapply(paraValues = paraValuesList, varPar = SigmaValuesList, FUN = function(paraValues, varPar) {
    if (!is.matrix(varPar)) {
      varPar <- varPar*diag(length(paraValues))
    }
    mvtnorm::dmvnorm(x = paraValues, sigma = varPar, log = TRUE)
  }))
}

setupVrecursionStep <- function(grid, v, baseVec1, baseVec2, K) {
  function(spaceTime1, spaceTime2) {
    if(!sameGridSection(spaceTime1, spaceTime2, grid)) {
      return(0)
    }
    V - baseVec1%*%K%*%baseVec2
  }
}

initVrecursion <- function(knots, covFct) {
  function(spaceTime1, spaceTime2) {
    initBList <- list(covFct(spaceTime1, knots), covFct(spaceTime2, knots))
    initInvK <- covFct(knots, knots)
    initK <- solve(initInvK)
    initV <- covFct(spaceTime1, spaceTime2)
    list(v = initV, K = initK, bList = initBList)
  }
}

setupFullrecursion <- function(gridList, knotsList, covFct) {
  initialVfun <- initVrecursion(knots = knotsList[[1]], covFct = covFct)
  function(spaceTime1, spaceTime2) {
    currentValues <- initialVfun(spaceTime1, spaceTime2)
    lapply(seq_along(gridList), FUN = function(resIndex) {
      incrementedVfun <- setupVrecursionStep(grid = gridList[[resIndex]], v = currentValues$v, baseVec1 = currentValues$bList[[1]], baseVec2 = currentValue$bList[[2]], K = currentValue$K)
      newB <- list(incrementedVfun(spaceTime1, knotsList))
    })
  }
}

buildRectanglePolygon <- function(corner1, corner2) {
  coords <- matrix(c(corner1[1], corner1[1], corner2[1], corner2[1], corner1[2], corner2[2], corner2[2], corner1[2]), nrow = 4, ncol = 2)
  Polygon(coors, hole = FALSE)
}

buildGrid <- function(spatialPointsGrid) {
  gridCoordAsList <- data.frame(t(spatialPointsGrid@data))
  rowIndices <- sort(unique(spatialPointsGrid@data[ , 1]))

  horizontalFunction <- function(horizontalIndex) {
    rowIndex <- rowIndices[horizontalIndex]
    subData <- spatialPointsGrid@data[spatialPointsGrid@data[ , 1] == rowIndex, ]
    colIndices <- sort(subData[ , 2])

    verticalFunction <- function(verticalIndex) {
      colIndex <- colIndices[verticalIndex]
      topLeftIndex <- match(data.frame(c(rowIndex, colIndex)), gridCoordAsList)
      topLeftCoord <- spatialPointsGrid@coords[topLeftIndex, ]
      topRight <- c(rowIndex, colIndices[verticalIndex+1])
      topRightIndex <- match(data.frame(topRight), gridCoordAsList)
      topRightCoord <- spatialPointsGrid@coords[topRightIndex, ]

      rowIndicesInColumn <- spatialPointsGrid@data[spatialPointsGrid@data[, 2] == topRight[[2]], 1]
      rowIndices
      bottomHorizontalIndex <-  max(match(rowIndex, sort(rowIndicesInColumn)) + 1)
      bottomRowIndex <- rowIndicesInColumn[[bottomHorizontalIndex]]
      bottomRight <- c(bottomRowIndex, topRight[[2]])
      bottomRightIndex <- match(data.frame(bottomRight), gridCoordAsList)
      bottomRightCoord <- spatialPointsGrid@coords[bottomRightIndex, ]
      bottomLeft <- c(bottomRowIndex, colIndex)
      bottomLeftIndex <- match(data.frame(bottomLeft), gridCoordAsList)
      bottomLeftCoord <- spatialPointsGrid@coords[bottomLeftIndex, ]
      polygonMatrix <- rbind(topLeftCoord, topRightCoord, bottomRightCoord, bottomLeftCoord, topLeftCoord)
      Polygon(polygonMatrix, hole = FALSE)
    }
    polygonsInRow <- lapply(1:(length(colIndices)-1), FUN = verticalFunction)
    polygonsInRow
  }

  allPolygons <- lapply(1:(length(rowIndices) - 1), FUN = horizontalFunction)
  SpatialPolygons(list(Polygons(do.call("c", allPolygons), ID = "ID")))
}


