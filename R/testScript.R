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

# Testing matrix inversion schemes

library(MASS)
library(Matrix)

k   <- 1000
rho <- .3

S       <- matrix(rep(rho, k*k), nrow=k)
diag(S) <- 1

dat <- mvrnorm(10000, mu=rep(0,k), Sigma=S) ### be patient!

R <- cor(dat)
Rsymm <- sparseMatrix(i= row(diag(k))[lower.tri(diag(k))], j= col(diag(k))[lower.tri(diag(k))], x = R[lower.tri(R)], symmetric=TRUE)
diag(Rsymm) <- 1

system.time(RI1 <- solve(R))
system.time(RI2 <- Matrix::chol2inv(Matrix::chol(R)))
system.time(RI3 <- qr.solve(R))
system.time(RI4 <- Matrix::chol2inv(Matrix::chol(Rsymm))) # Slower than RI2.

all.equal(RI1, RI2)
all.equal(RI1, RI3)

