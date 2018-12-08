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

List logLikMRAcpp(SEXP treePointer, NumericVector & covParameters, NumericVector & fixedEffectParameters, double & errorSD, double & fixedEffSD)
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
    pointedTree->CenterResponse() ;

    pointedTree->ComputeMRAloglik() ;
    // pointedTree->ComputeGlobalLogLik() ;
    logLikVal = pointedTree->GetMRAlogLik() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return List::create(Named("logLik") = logLikVal, Named("gridPointer") = treePointer) ;
}

// [[Rcpp::export]]

List predictMRArcpp(SEXP treePointer, NumericMatrix predSpatialCoor, IntegerVector predTime, NumericMatrix covariatesAtPredLocs) {
  uvec predictionTime = as<uvec>(predTime) ;
  mat predSp = as<mat>(predSpatialCoor) ;
  mat predCovariates = as<mat>(covariatesAtPredLocs) ;
  std::vector<GaussDistParas> predsInZones ;
  std::vector<vec> meanVecVec ;
  std::vector<mat> covMatVec ;
  vec predResponses = zeros<vec>(predictionTime.size()) ;
  inputdata dataForPreds(predResponses, predSp, predictionTime, predCovariates) ;

  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    spatialcoor predLocs(predSp, predictionTime) ;
    predsInZones = pointedTree->ComputeConditionalPrediction(dataForPreds) ;
    cout << "Entering copy loop... \n" ;
    printf("Number of zones: %u \n", predsInZones.size()) ;
    cout << "Indexed value: " << predsInZones.at(0).covPara(0,0) << "\n" ;
    uint counter = 0 ;
    for (auto & i : predsInZones) {
      if (i.covPara(0,0) != -1) { // A value under 0 on the diagonal indicates that the region had no observations for prediction
        meanVecVec.push_back(i.meanPara) ;
        covMatVec.push_back(i.covPara) ;
        printf("Counter value: %u \n", counter) ;
        counter += 1 ;
      }
    }
    cout << "Leaving copy loop... \n" ;
    pointedTree->CleanPredictionComponents() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return List::create(Named("means") = Rcpp::wrap(meanVecVec),
                      Named("covMatricesNoDim") = Rcpp::wrap(covMatVec),
                      Named("gridPointer") = treePointer) ;
}

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
