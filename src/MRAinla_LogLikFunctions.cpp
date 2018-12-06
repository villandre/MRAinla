// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
// #include "gperftools/profiler.h"

#include "AugTree.h"

using namespace arma;
using namespace Rcpp;
using namespace MRAinla;

typedef unsigned int uint ;

// // [[Rcpp::export]]
// SEXP start_profiler(SEXP str) {
//   ProfilerStart(as<const char*>(str));
//   return R_NilValue;
// }

// // [[Rcpp::export]]
// SEXP stop_profiler() {
//   ProfilerStop();
//   return R_NilValue;
// }

// [[Rcpp::export]]

SEXP setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, IntegerVector obsTime,
                  NumericMatrix covariateMatrix, uint M, NumericVector lonRange, NumericVector latRange,
                  IntegerVector timeRange, uint randomSeed, uint cutForTimeSplit)
{
  vec lonR = as<vec>(lonRange) ;
  vec latR = as<vec>(latRange) ;
  uvec timeR = as<uvec>(timeRange) ;
  vec response = as<vec>(responseValues) ;
  mat sp = as<mat>(spCoords) ;
  uvec time = as<uvec>(obsTime) ;
  unsigned long int seedForRNG = randomSeed ;
  mat covariateMat = as<mat>(covariateMatrix) ;

  AugTree * MRAgrid = new AugTree(M, lonR, latR, timeR, response, sp, time, cutForTimeSplit, seedForRNG, covariateMat) ;

  XPtr<AugTree> p(MRAgrid, false) ; // Disabled automatic garbage collection.

  return List::create(Named("gridPointer") = p) ;
}

// [[Rcpp::export]]

List logLikCpp(SEXP treePointer, NumericVector & covParameters, NumericVector & fixedEffectParameters, double & errorSD, double & fixedEffSD)
{
  //omp_set_num_threads(numOpenMP) ;
  double logLikVal = 0;
  vec fixedEffVec = as<vec>(fixedEffectParameters) ;
  vec covParVec = as<vec>(covParameters) ;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    pointedTree->SetCovParameters(covParVec) ;
    pointedTree->SetFixedEffParameters(fixedEffVec) ;
    pointedTree->SetErrorSD(errorSD) ;
    pointedTree->SetFixedEffSD(fixedEffSD) ;

    pointedTree->ComputeMRAloglik() ;
    pointedTree->ComputeGlobalLogLik() ;
    logLikVal = pointedTree->GetMRAlogLik() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return List::create(Named("logLik") = logLikVal, Named("gridPointer") = treePointer) ;
}

// [[Rcpp::export]]

List testFunction() {
  std::vector<mat> myVec(10) ;
  for(auto & i : myVec) {
    i.resize(2,4) ;
    i.fill(0) ;
  }
  List returnList = Rcpp::wrap(myVec) ;
  return returnList ;
}


// List predictMRA(SEXP treePointer, NumericMatrix predSpatialCoor, IntegerVector predTime) {
//   uvec predictionTime = as<uvec>(predTime) ;
//   mat predSp = as<mat>(predSpatialCoor) ;
//   std::vector<GaussDistParas> predsInZones ;
//
//   if (!(treePointer == NULL))
//   {
//     XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
//     spatialcoor predLocs(predSp, predictionTime) ;
//     predsInZones = pointedTree->ComputeConditionalPrediction(predLocs) ;
//   }
//   else
//   {
//     throw Rcpp::exception("Pointer to MRA grid is null." ) ;
//   }
//   List convert(predInZones)
// }

// [[Rcpp::export]]

NumericMatrix inla(SEXP treePointer, NumericMatrix predictionLocations, IntegerVector predictionTime, double stepSize) {

  arma::mat posteriorMatrix ;

  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    spatialcoor predictionST(as<mat>(predictionLocations), as<uvec>(predictionTime)) ;
    posteriorMatrix = pointedTree->ComputePosteriors(predictionST, stepSize) ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return Rcpp::wrap(posteriorMatrix) ;
}
