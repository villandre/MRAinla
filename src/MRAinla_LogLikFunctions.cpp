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

  return List::create(Named("gridPointer") = p) ;
}

// [[Rcpp::export]]

List logLikCpp(SEXP treePointer)
{
  //omp_set_num_threads(numOpenMP) ;
  double logLikVal = 0;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    pointedTree->ComputeLoglik() ;
    double logLikVal = pointedTree->GetLoglik() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return List::create(Named("logLik") = logLikVal, Named("gridPointer") = treePointer) ;
}
