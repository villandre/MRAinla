// [[Rcpp::depends(RcppArmadillo)]]

#include <iostream>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>
#include "gperftools/profiler.h"

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
                  NumericVector timeRange, uint randomSeed, uint cutForTimeSplit)
{
  fvec lonR = as<fvec>(lonRange) ;
  fvec latR = as<fvec>(latRange) ;
  fvec timeR = as<fvec>(timeRange) ;
  vec response = as<vec>(responseValues) ;
  fmat sp = as<fmat>(spCoords) ;
  fvec time = as<fvec>(obsTime) ;
  unsigned long int seedForRNG = randomSeed ;
  fmat covariateMat = as<fmat>(covariateMatrix) ;

  AugTree * MRAgrid = new AugTree(M, lonR, latR, timeR, response, sp, time, cutForTimeSplit, seedForRNG, covariateMat) ;

  XPtr<AugTree> p(MRAgrid, false) ; // Disabled automatic garbage collection.

  return List::create(Named("gridPointer") = p) ;
}

// // [[Rcpp::export]]
//
// List logLikMRAcpp(SEXP treePointer, NumericVector & covParameters, NumericVector & fixedEffectParameters,
//                   double & errorSD, double & fixedEffSD, NumericVector & fieldValues)
// {
//   //omp_set_num_threads(numOpenMP) ;
//   double logLikVal = 0;
//   vec fixedEffVec = as<vec>(fixedEffectParameters) ;
//   vec covParVec = as<vec>(covParameters) ;
//   vec fieldValuesVec = as<vec>(fieldValues) ;
//   if (!(treePointer == NULL))
//   {
//     XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
//     bool testEqual = false ;
//     if (covParVec.size() == pointedTree->GetCovParameters().size()) {
//       testEqual = approx_equal(covParVec, pointedTree->GetCovParameters(), "absdiff", 1e-6) ;
//     } // If covariance parameters haven't changed since the last evaluation, there's no need to re-run functions involving only covariance calculations.
//
//     pointedTree->SetCovParameters(covParVec) ;
//     pointedTree->SetFixedEffParameters(fixedEffVec) ;
//     pointedTree->SetErrorSD(errorSD) ;
//     pointedTree->SetFixedEffSD(fixedEffSD) ;
//     // pointedTree->CenterResponse() ;
//
//     logLikVal = pointedTree->ComputeMRAloglik(fieldValuesVec, testEqual) ;
//     // pointedTree->ComputeGlobalLogLik() ;
//   }
//   else
//   {
//     throw Rcpp::exception("Pointer to MRA grid is null." ) ;
//   }
//   return List::create(Named("logLik") = logLikVal, Named("gridPointer") = treePointer) ;
// }

// // [[Rcpp::export]]

// List predictMRArcpp(SEXP treePointer, NumericMatrix predSpatialCoor, IntegerVector predTime, NumericMatrix covariatesAtPredLocs) {
//   uvec predictionTime = as<uvec>(predTime) ;
//   mat predSp = as<mat>(predSpatialCoor) ;
//   mat predCovariates = as<mat>(covariatesAtPredLocs) ;
//   std::vector<GaussDistParas> predsInZones ;
//   std::vector<vec> meanVecVec ;
//   std::vector<mat> covMatVec ;
//   vec predResponses = zeros<vec>(predictionTime.size()) ;
//   inputdata dataForPreds(predResponses, predSp, predictionTime, predCovariates) ;
//
//   if (!(treePointer == NULL))
//   {
//     XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
//     spatialcoor predLocs(predSp, predictionTime) ;
//     predsInZones = pointedTree->ComputeConditionalPrediction(dataForPreds) ;
//
//     for (auto & i : predsInZones) {
//       if (i.covPara(0,0) != -1) { // A value under 0 on the diagonal indicates that the region had no observations for prediction
//         meanVecVec.push_back(i.meanPara) ;
//         covMatVec.push_back(i.covPara) ;
//       }
//     }
//     pointedTree->CleanPredictionComponents() ;
//   }
//   else
//   {
//     throw Rcpp::exception("Pointer to MRA grid is null." ) ;
//   }
//   return List::create(Named("means") = Rcpp::wrap(meanVecVec),
//                       Named("covMatricesNoDim") = Rcpp::wrap(covMatVec),
//                       Named("gridPointer") = treePointer) ;
// }

// // [[Rcpp::export]]
//
// NumericMatrix inla(SEXP treePointer, NumericMatrix predictionLocations, IntegerVector predictionTime, double stepSize) {
//
//   arma::mat posteriorMatrix ;
//
//   if (!(treePointer == NULL))
//   {
//     XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
//     spatialcoor predictionST(as<mat>(predictionLocations), as<uvec>(predictionTime)) ;
//     posteriorMatrix = pointedTree->ComputePosteriors(predictionST, stepSize) ;
//   }
//   else
//   {
//     throw Rcpp::exception("Pointer to MRA grid is null." ) ;
//   }
//   return Rcpp::wrap(posteriorMatrix) ;
// }

// [[Rcpp::export]]

double funForOptimJointHyperMarginal(SEXP treePointer, Rcpp::NumericVector MRAhyperparas,
         double fixedEffSD, double errorSD, Rcpp::List MRAcovParasIGalphaBeta,
         NumericVector fixedEffIGalphaBeta, NumericVector errorIGalphaBeta) {
  arma::mat posteriorMatrix ;
  double outputValue = 0 ;

  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    // outputValue = pointedTree->ComputeJointPsiMarginalPropConstant(as<vec>(MRAhyperparas),
    //                                              fixedEffSD, errorSD, hyperAlpha, hyperBeta) ;
    // ProfilerStart("/home/luc/Downloads/myprofile.log") ;
    // The alpha's and beta's for the IG distribution of the hyperparameters do not change.
    if (pointedTree->GetMRAcovParasIGalphaBeta().size() == 0) {
      std::vector<IGhyperParas> MRAalphaBeta ;
      for (auto & i : MRAcovParasIGalphaBeta) {
        vec alphaBetaVec = Rcpp::as<vec>(i) ;
        IGhyperParas alphaBetaStruct(alphaBetaVec.at(0), alphaBetaVec.at(1)) ;
        MRAalphaBeta.push_back(alphaBetaStruct) ;
      }
      pointedTree->SetMRAcovParasIGalphaBeta(MRAalphaBeta) ;
      vec fixedEffAlphaBeta = Rcpp::as<vec>(fixedEffIGalphaBeta) ;
      vec errorAlphaBeta = Rcpp::as<vec>(errorIGalphaBeta) ;

      pointedTree->SetFixedEffIGalphaBeta(IGhyperParas(fixedEffAlphaBeta.at(0), fixedEffAlphaBeta.at(1))) ;
      pointedTree->SetErrorIGalphaBeta(IGhyperParas(errorAlphaBeta.at(0), errorAlphaBeta.at(1))) ;
    }
    pointedTree->SetErrorSD(errorSD) ;
    pointedTree->SetFixedEffSD(fixedEffSD) ;
    pointedTree->SetMRAcovParas(MRAhyperparas) ;

    outputValue = pointedTree->ComputeJointPsiMarginal() ;
    // ProfilerStop() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  printf("Marginal joint Psi: %.4e \n", outputValue) ;
  return outputValue ;
}

