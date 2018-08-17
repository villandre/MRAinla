library(sp)

pointsIndices <- cbind(c(1,1,2,2,2,3,3,3),c(1,3,1,2,3,1,2,3))
pointsCoords <- pointsIndices/5

spatialGrid <- SpatialPointsDataFrame(coords = pointsCoords, data = as.data.frame(pointsIndices))

buildGrid(spatialGrid)
