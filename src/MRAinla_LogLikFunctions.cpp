// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

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

List setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, NumericVector obsTime,
                  NumericMatrix covariateMatrix, uint M, NumericVector lonRange, NumericVector latRange,
                  NumericVector timeRange, uint randomSeed, uint cutForTimeSplit, bool splitTime)
{
  vec lonR = as<vec>(lonRange) ;
  vec latR = as<vec>(latRange) ;
  vec timeR = as<vec>(timeRange) ;
  vec response = as<vec>(responseValues) ;
  mat sp = as<mat>(spCoords) ;
  vec time = as<vec>(obsTime) ;

  unsigned long int seedForRNG = randomSeed ;

  mat covariateMat = as<mat>(covariateMatrix) ;

  AugTree * MRAgrid = new AugTree(M, lonR, latR, timeR, response, sp, time, cutForTimeSplit, seedForRNG, covariateMat, splitTime) ;

  XPtr<AugTree> p(MRAgrid, false) ; // Disabled automatic garbage collection.

  return List::create(Named("gridPointer") = p) ;
}

// [[Rcpp::export]]

double LogJointHyperMarginalToWrap(SEXP treePointer, Rcpp::List MRAhyperparas,
         double fixedEffSD, double errorSD, Rcpp::List MRAcovParasGammaAlphaBeta,
         Rcpp::NumericVector FEmuVec, NumericVector fixedEffGammaAlphaBeta,
         NumericVector errorGammaAlphaBeta, bool matern, double spaceNuggetSD, double timeNuggetSD,
         bool recordFullConditional, Rcpp::Function optimFun, Rcpp::Function LUfun) {
  arma::mat posteriorMatrix ;
  double outputValue = 0 ;

  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.

    // The alpha's and beta's for the gamma distribution of the hyperparameters do not change.
    if (!pointedTree->CheckMRAcovParasGammaAlphaBeta()) {
      pointedTree->ToggleGammaParasSet() ;
      pointedTree->SetMRAcovParasGammaAlphaBeta(MRAcovParasGammaAlphaBeta) ;

      vec fixedEffAlphaBeta = Rcpp::as<vec>(fixedEffGammaAlphaBeta) ;
      vec errorAlphaBeta = Rcpp::as<vec>(errorGammaAlphaBeta) ;
      vec FEmu = Rcpp::as<vec>(FEmuVec) ;

      pointedTree->SetFixedEffGammaAlphaBeta(GammaHyperParas(fixedEffAlphaBeta.at(0), fixedEffAlphaBeta.at(1))) ;
      pointedTree->SetErrorGammaAlphaBeta(GammaHyperParas(errorAlphaBeta.at(0), errorAlphaBeta.at(1))) ;
      pointedTree->SetFEmu(FEmu) ;
      pointedTree->SetMatern(matern) ;
      pointedTree->SetSpaceAndTimeNuggetSD(spaceNuggetSD, timeNuggetSD) ;
      std::vector<TreeNode *> tipNodes = pointedTree->GetLevelNodes(pointedTree->GetM()) ;
      for (auto & i : tipNodes) {
        i->SetUncorrSD(0.05) ;
      }
    }
    pointedTree->SetErrorSD(errorSD) ;
    pointedTree->SetFixedEffSD(fixedEffSD) ;

    pointedTree->SetMRAcovParas(MRAhyperparas) ;
    pointedTree->SetRecordFullConditional(recordFullConditional) ;
    // ProfilerStart("/home/luc/Downloads/myprofile.log") ;

    pointedTree->ComputeLogJointPsiMarginal(optimFun, LUfun) ;

    // ProfilerStop() ;
    // throw Rcpp::exception("Stop for profiling... \n") ;
    outputValue = pointedTree->GetLogJointPsiMarginal() ;
    // printf("Marginal joint Psi: %.4e \n \n \n", pointedTree->GetLogJointPsiMarginal()) ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return outputValue ;
}

// [[Rcpp::export]]

arma::vec GetFullCondMean(SEXP treePointer) {
  vec outputVec ;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    outputVec = pointedTree->GetFullCondMean() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return outputVec ;
}

// [[Rcpp::export]]

arma::vec GetFullCondSDs(SEXP treePointer) {
  vec outputVec ;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    outputVec = pointedTree->GetFullCondSDs() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return outputVec ;
}

// [[Rcpp::export]]

Rcpp::List ComputeCondPredStats(SEXP treePointer, NumericMatrix spCoordsForPredict, NumericVector timeForPredict,
                                 NumericMatrix covariateMatrixForPredict) {
  sp_mat Hmat, Hmean, HmeanSq ;
  mat HmeanMat ;
  vec Evar ;

  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    mat spCoords = Rcpp::as<mat>(spCoordsForPredict) ;
    vec time = Rcpp::as<vec>(timeForPredict) ;
    mat covariates = Rcpp::as<mat>(covariateMatrixForPredict) ;
    // ProfilerStart("/home/luc/Downloads/myprofile.log") ;

    Hmat = conv_to<mat>::from(pointedTree->ComputeHpred(spCoords, time, covariates)) ;
    sp_mat Hmean = Hmat * conv_to<sp_mat>::from(pointedTree->GetFullCondMean()) ;
    HmeanMat = conv_to<mat>::from(Hmean) ;
    Evar = pointedTree->ComputeEvar(Hmat) ;

    // ProfilerStop() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return Rcpp::List::create(Named("Hmean") = HmeanMat, Named("Evar") = Evar) ;
}
