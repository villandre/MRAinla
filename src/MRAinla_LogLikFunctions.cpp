// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

typedef unsigned int uint ;

// // [[Rcpp::export]]
//
// SEXP setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, DateVector obsTime, uint M, NumericVector lonRange, NumericVector latRange, DateVector timeRange, uint randomSeed, uint cutForTimeSplit = 200, uint cutForLatSplit = 200)
// {
//   // TO_DO
// }
//
// // [[Rcpp::export]]
//
// List logLikCpp(SEXP treePointer, uint brickDepth)
// {
//   // TO_DO
// }
