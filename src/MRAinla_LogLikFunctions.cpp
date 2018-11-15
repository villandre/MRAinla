// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>

#include "AugTree.h"

using namespace arma;
using namespace Rcpp;
using namespace MRAinla;

typedef unsigned int uint ;

// [[Rcpp::export]]

SEXP setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, IntegerVector obsTime,
                  uint M, NumericVector lonRange, NumericVector latRange,
                  IntegerVector timeRange, uint randomSeed, uint cutForTimeSplit,
                  NumericVector covarianceParameter)
{
  vec lonR = as<vec>(lonRange) ;
  vec latR = as<vec>(latRange) ;
  uvec timeR = as<uvec>(timeRange) ;
  vec response = as<vec>(responseValues) ;
  mat sp = as<mat>(spCoords) ;
  uvec time = as<uvec>(obsTime) ;
  vec covParaVec = as<vec>(covarianceParameter) ;
  unsigned long int seedForRNG = randomSeed ;

  AugTree * MRAgrid = new AugTree(M, lonR, latR, timeR, response, sp, time, cutForTimeSplit, seedForRNG, covParaVec) ;

  XPtr<AugTree> p(MRAgrid, false) ; // Disabled automatic garbage collection.

  return List::create(Named("logLik") = 1,
                      Named("gridPointer") = p) ;
}
//
// // [[Rcpp::export]]
//
// List logLikCpp(SEXP treePointer, uint brickDepth)
// {
//   // TO_DO
// }
