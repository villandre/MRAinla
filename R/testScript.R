library(sp)

pointsIndices <- cbind(c(1,1,2,2,2,3,3,3),c(1,3,1,2,3,1,2,3))
pointsCoords <- pointsIndices/5

spatialGrid <- SpatialPointsDataFrame(coords = pointsCoords, data = as.data.frame(pointsIndices))

foo <- buildGrid(spatialGrid)

## Testing the grid functions

aGrid <- gridConstructor(lonBreaks = 1:9, latBreaks = 1:5, timeBreaks = 1:19, lonExtent = c(0,10), latExtent = c(0,6), timeExtent = c(0,20))
aGrid # print method is ok.

set.seed(10)
sampleObservations <- SpatialPoints(cbind(runif(50)*10, runif(50)*6))
sampleKnots <- SpatialPoints(cbind(runif(5)*10, runif(5)*6))

plot(aGrid, observationsAsPoints = sampleObservations, knotsAsPoints = sampleKnots)

