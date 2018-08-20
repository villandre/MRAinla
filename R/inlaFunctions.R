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
    v - baseVec1%*%K%*%baseVec2
  }
}

sameGridSection <- function(spaceTime1, spaceTime2, grid, timeBreaks) {
  timeOrder <- order(c(spaceTime1$time, spaceTime2$time, timeBreaks))
  if (abs(timeOrder[1]-timeOrder[2]) > 1) {
    return(FALSE)
  }
  polyForPoints <- sapply(grid@polygons, FUN = function(aPolygon) {
    pointsInPolyTest <- point.in.polygon(point.x = c(spaceTime1$coords[[1]], spaceTime2$coords[[1]]) , point.y = c(spaceTime1$coords[[2]], spaceTime2$coords[[2]]), pol.x = aPolygon@coords[ , 1], pol.y = aPolygon@coords[ , 2])
    pointsInPolyTest <- replace(pointsInPolyTest, pointsInPolyTest > 1, 1)
    if (pointsInPolyTest[[1]] != pointsInPolyTest[[2]]) {
      return(FALSE)
    }
    TRUE
  })
  if (!all(polyForPoints)) {
    return(FALSE)
  }
  TRUE
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
      newB <- list(incrementedVfun(spaceTime1, knotsList), incrementedVfun(spaceTime2, knotsList))
      newV <- incrementedVfun(spaceTime1, spaceTime2)

    })
  }
}

buildRectanglePolygon <- function(corner1, corner2) {
  coords <- matrix(c(corner1[1], corner1[1], corner2[1], corner2[1], corner1[2], corner2[2], corner2[2], corner1[2]), nrow = 4, ncol = 2)
  Polygon(coors, hole = FALSE)
}

# Alternatively, the grids could be represented as a RasterBrick object.

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

getPointsFromRaster2D <- function(rasterMap, rasterValue, nodeCoords) {
  cellsInRegion <- which(raster::getValues(rasterMap) == rasterValue)
  extentObject <- raster::extentFromCells(rasterMap, cells = cellsInRegion)
  raster::intersect(nodeCoords, extentObject)
}

getPointsFromRasterSpatiotemp <- function(rasterBrick, rasterValuesVec, nodeCoordsList) {
  lapply(seq_along(rasterValuesVec), FUN = function(timeIndex) {
    getPointsFromRaster2D(subset(rasterBrick, timeIndex), rasterValue = rasterValuesVec[[timeIndex]], nodeCoords =  nodeCoordsList[[timeIndex]])
  })
}

mergeGridSections <- function(rasterMap, sectionIndex1, sectionIndex2) {
  values(rasterMap) <- replace(getValues(rasterMap), which(getValues(rasterMap) == max(c(sectionIndex1, sectionIndex2))), min(c(sectionIndex1, sectionIndex2)))
  values(rasterMap)[which(getValues(rasterMap) > max(c(sectionIndex1, sectionIndex2))))] <- getValues(rasterMap)[which(getValues(rasterMap) > max(c(sectionIndex1, sectionIndex2))))] - 1
  rasterMap
}

# The raster bricks also have an associated matrix with M columns indicating the layer breakdown at each resolution, e.g. for M = 2, we have this matrix equal to cbind(c(1,1,1,1), c(1,1,2,2), c(1,2,3,4)) (for resolutions 0,1, and 2).

# A 3D grid can be represented parsemoniously by a list of 3 vectors, one for each dimension. We assume length, width, and depth are all split into small pixels. Pixels can also be rectangles. A split translates into a change into one of the vectors. The object needs a range component to help associate pixels with spatiotemporal coordinates.

gridConstructor <- function(lonBreaks, latBreaks, timeBreaks, lonExtent, latExtent, timeExtent) {
  grids <- mapply(breaks = list(lonBreaks, latBreaks, timeBreaks), extent = list(lonExtent, latExtent, timeExtent), FUN = function(breaks, extent) {
    list(breaks = breaks, extent = extent)
  }, SIMPLIFY = FALSE)
  names(grids) <- c("longitude", "latitude", "time")
  class(grids) <- "Spacetimegrid"
  grids
}

splitGridSection <- function(gridObject, dimension = c("lon", "lat", "time"), breakPos) {
  gridObject[[dimension]]$breaks <- sort(c(gridObject[[dimension]]$breaks, breakPos))
  gridObject
}

print.Spacetimegrid <- function(x) {
  cat("Spacetime grid object \n")
  cat("Longitude extent: ", x$longitude$extent , "\n")
  cat("Latitude extent: ", x$latitude$extent, "\n")
  cat("Time extent: ", x$time$extent, "\n \n")
  cat("Grid dimensions: (", length(x$longitude$breaks) + 1, " ", length(x$latitude$breaks)+1, " ", length(x$time$breaks) + 1, ")\n")
  cat("Note: This gives the number of regions, which may have different volumes.\n\n")
  invisible(NULL)
}

# This function actually plots a slice from a Spacetimegrid object.

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
    combinedSP <- bind(combinedSP, knotSP)
    plotColours <- combinedSP@data$identity
  }
  plot(combinedSP, col = as.character(plotColours), xlab = "Longitude", ylab = "Latitude", xlim = x$longitude$extent, ylim = x$latitude$extent, cex = 0.8, pch = 18)
  abline(v = c(x$longitude$extent[1], x$longitude$breaks, x$longitude$extent[2]), lty = 3, col = "blue4")
  abline(h = c(x$latitude$extent[1], x$latitude$breaks, x$latitude$extent[2]), lty = 3, col = "blue4")
}



