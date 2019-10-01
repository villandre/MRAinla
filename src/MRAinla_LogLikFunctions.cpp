// [[Rcpp::depends(RcppEigen)]]

#include <iostream>
#include <string>
#include <algorithm>
#include <RcppEigen.h>
// #include "gperftools/profiler.h"

#include "AugTree.h"

using namespace Eigen;
using namespace Rcpp;
using namespace MRAinla;

// [[Rcpp::export]]

List setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, NumericMatrix predCoords,
                  NumericVector obsTime,  NumericVector predTime, NumericMatrix covariateMatrix,
                  NumericMatrix predCovariateMatrix, uint M,
                  NumericVector lonRange, NumericVector latRange, NumericVector timeRange,
                  uint randomSeed, uint cutForTimeSplit, bool splitTime,
                  int numKnotsRes0, int J)
{
  ArrayXd lonRinit = as<ArrayXd>(lonRange) ;
  ArrayXd latRinit = as<ArrayXd>(latRange) ;
  ArrayXd timeRinit = as<ArrayXd>(timeRange) ;
  Array2d lonR = lonRinit.segment(0,2) ;
  Array2d latR = latRinit.segment(0,2) ;
  Array2d timeR = timeRinit.segment(0,2) ;
  vec response = as<vec>(responseValues) ;
  ArrayXXd sp = as<ArrayXXd>(spCoords) ;
  ArrayXXd predSp = as<ArrayXXd>(predCoords) ;
  ArrayXd predTimeVec = as<ArrayXd>(predTime) ;
  ArrayXXd predCovariates = as<ArrayXXd>(predCovariateMatrix) ;
  ArrayXd time = as<ArrayXd>(obsTime) ;

  unsigned long int seedForRNG = randomSeed ;

  ArrayXXd covariateMat = as<ArrayXXd>(covariateMatrix) ;

  AugTree * MRAgrid = new AugTree(M, lonR, latR, timeR, response, sp, time, predCovariates, predSp, predTimeVec, cutForTimeSplit, seedForRNG, covariateMat, splitTime, numKnotsRes0, J) ;

  XPtr<AugTree> p(MRAgrid, false) ; // Disabled automatic garbage collection.

  return List::create(Named("gridPointer") = p) ;
}

// [[Rcpp::export]]

double LogJointHyperMarginalToWrap(SEXP treePointer, Rcpp::List MRAhyperparas,
         double fixedEffSD, double errorSD, Rcpp::List MRAcovParasGammaAlphaBeta,
         Rcpp::NumericVector FEmuVec, NumericVector fixedEffGammaAlphaBeta,
         NumericVector errorGammaAlphaBeta, bool matern, double spaceNuggetSD, double timeNuggetSD,
         bool recordFullConditional) {
  mat posteriorMatrix ;
  double outputValue = 0 ;

  // if (!(treePointer == NULL))
  // {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.

    // The alpha's and beta's for the gamma distribution of the hyperparameters do not change.
    if (!pointedTree->CheckMRAcovParasGammaAlphaBeta()) {
      pointedTree->ToggleGammaParasSet() ;
      pointedTree->SetMRAcovParasGammaAlphaBeta(MRAcovParasGammaAlphaBeta) ;

      vec fixedEffAlphaBeta = Rcpp::as<vec>(fixedEffGammaAlphaBeta) ;
      vec errorAlphaBeta = Rcpp::as<vec>(errorGammaAlphaBeta) ;
      vec FEmu = Rcpp::as<vec>(FEmuVec) ;

      pointedTree->SetFixedEffGammaAlphaBeta(GammaHyperParas(fixedEffAlphaBeta(0), fixedEffAlphaBeta(1))) ;
      pointedTree->SetErrorGammaAlphaBeta(GammaHyperParas(errorAlphaBeta(0), errorAlphaBeta(1))) ;
      pointedTree->SetFEmu(FEmu) ;
      pointedTree->SetMatern(matern) ;
      pointedTree->SetSpaceAndTimeNuggetSD(spaceNuggetSD, timeNuggetSD) ;
      std::vector<TreeNode *> tipNodes = pointedTree->GetLevelNodes(pointedTree->GetM()) ;
      for (auto & i : tipNodes) {
        i->SetUncorrSD(0.001) ; // Is this a nugget effect?
      }
    }
    pointedTree->SetErrorSD(errorSD) ;
    pointedTree->SetFixedEffSD(fixedEffSD) ;

    pointedTree->SetMRAcovParas(MRAhyperparas) ;
    pointedTree->SetRecordFullConditional(recordFullConditional) ;

    pointedTree->ComputeLogJointPsiMarginal() ;

    outputValue = pointedTree->GetLogJointPsiMarginal() ;
    // printf("Marginal joint Psi: %.4e \n \n \n", pointedTree->GetLogJointPsiMarginal()) ;
  // }
  // else
  // {
  //   throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  // }
  return outputValue ;
}

// [[Rcpp::export]]

Eigen::VectorXd GetFullCondMean(SEXP treePointer) {
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

Eigen::VectorXd GetFullCondSDs(SEXP treePointer) {
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
                                 NumericMatrix covariateMatrixForPredict, int batchSize) {
  sp_mat Hmat, Hmean, HmeanSq ;
  mat HmeanMat ;
  vec Evar ;
  std::cout << "Entered ComputeCondPredStats \n" ;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
    mat spCoords = Rcpp::as<mat>(spCoordsForPredict) ;
    vec time = Rcpp::as<vec>(timeForPredict) ;
    mat covariates = Rcpp::as<mat>(covariateMatrixForPredict) ;

    vec Hmean = pointedTree->GetHmatPred() * pointedTree->GetFullCondMean() ;
    std::cout << "Computing Evar! \n" ;
    Evar = pointedTree->ComputeEvar(batchSize) ;
    std::cout << "Done! \n" ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return Rcpp::List::create(Named("Hmean") = Hmean, Named("Evar") = Evar) ;
}

// [[Rcpp::export]]

int GetNumTips(SEXP treePointer) {
  XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
  int value = pointedTree->GetNumTips() ;
  return value ;
}

// [[Rcpp::export]]

Eigen::VectorXi GetPredObsOrder(SEXP treePointer) {
  XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
  uvec value = pointedTree->GetObsOrderForHpredMat() ;
  return value ;
}
