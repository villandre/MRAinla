#include <iostream>
#include <string>
#include <algorithm>
#include <RcppEigen.h>
// #include "gperftools/profiler.h"

#include "Forest.h"

using namespace Eigen;
using namespace Rcpp;
using namespace MRAinla;
using namespace std;

// [[Rcpp::export]]

List setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, NumericMatrix predCoords,
                  NumericVector obsTime,  NumericVector predTime, NumericMatrix covariateMatrix,
                  NumericMatrix predCovariateMatrix, uint Mlon, uint Mlat,
                  NumericVector lonRange, NumericVector latRange,
                  uint randomSeed, int numKnotsRes0, double J, String distMethod)
{
  ArrayXd lonRinit = as<ArrayXd>(lonRange) ;
  ArrayXd latRinit = as<ArrayXd>(latRange) ;
  Array2d lonR = lonRinit.segment(0,2) ;
  Array2d latR = latRinit.segment(0,2) ;
  vec response = as<vec>(responseValues) ;
  ArrayXd responseArray = as<Eigen::ArrayXd>(responseValues) ;
  ArrayXXd sp = as<ArrayXXd>(spCoords) ;
  ArrayXXd predSp = as<ArrayXXd>(predCoords) ;
  ArrayXd predTimeVec = as<ArrayXd>(predTime) ;
  ArrayXXd predCovariates = as<ArrayXXd>(predCovariateMatrix) ;
  ArrayXd time = as<ArrayXd>(obsTime) ;
  string dMethod = as<std::string>(Rcpp::wrap(distMethod)) ;

  unsigned long int seedForRNG = randomSeed ;

  ArrayXXd covariateMat = as<ArrayXXd>(covariateMatrix) ;

  Forest * MRAgrids = new Forest(Mlon, Mlat, lonR, latR, response, sp, time, predCovariates, predSp, predTimeVec, seedForRNG, covariateMat, numKnotsRes0, J, dMethod) ;
  XPtr<Forest> p(MRAgrids, false) ; // Disabled automatic garbage collection.

  return List::create(Named("gridPointer") = p) ;
}

// [[Rcpp::export]]

double LogJointHyperMarginalToWrap(SEXP forestPointer, Rcpp::List MRAhyperparas,
         double timeCovPara, double fixedEffSD, double errorSD, Rcpp::List MRAcovParasGammaAlphaBeta,
         Rcpp::NumericVector FEmuVec, NumericVector fixedEffGammaAlphaBeta,
         NumericVector errorGammaAlphaBeta, NumericVector timeGammaAlphaBeta,
         double spaceNuggetSD, bool recordFullConditional, bool processPredictions) {
  mat posteriorMatrix ;
  double outputValue = 0 ;

  // if (!(forestPointer == NULL))
  // {
    XPtr<Forest> pointedForest(forestPointer) ; // Becomes a regular pointer again.

    // The alpha's and beta's for the gamma distribution of the hyperparameters do not change.

    if (!pointedForest->CheckMRAcovParasGammaAlphaBeta()) {
      pointedForest->ToggleGammaParasSet() ;
      pointedForest->SetMRAcovParasGammaAlphaBeta(MRAcovParasGammaAlphaBeta) ;

      vec fixedEffAlphaBeta = Rcpp::as<vec>(fixedEffGammaAlphaBeta) ;
      vec errorAlphaBeta = Rcpp::as<vec>(errorGammaAlphaBeta) ;
      vec FEmu = Rcpp::as<vec>(FEmuVec) ;

      pointedForest->SetFixedEffGammaAlphaBeta(GammaHyperParas(fixedEffAlphaBeta(0), fixedEffAlphaBeta(1))) ;
      pointedForest->SetErrorGammaAlphaBeta(GammaHyperParas(errorAlphaBeta(0), errorAlphaBeta(1))) ;
      pointedForest->SetTimeCovParaGammaAlphaBeta(GammaHyperParas(timeGammaAlphaBeta(0), timeGammaAlphaBeta(1))) ;
      pointedForest->SetFEmu(FEmu) ;
      pointedForest->SetSpaceNuggetSD(spaceNuggetSD) ;
      // std::vector<TreeNode *> tipNodes = pointedForest->GetLevelNodes(pointedForest->GetM()) ;
      // for (auto & i : tipNodes) {
      //   i->SetUncorrSD(0.001) ; // Is this a nugget effect?
      // }
    }
    pointedForest->SetErrorSD(errorSD) ;
    pointedForest->SetFixedEffSD(fixedEffSD) ;

    pointedForest->SetMRAcovParas(MRAhyperparas) ;
    pointedForest->SetRecordFullConditional(recordFullConditional) ;
    pointedForest->SetProcessPredictions(processPredictions) ;

    pointedForest->ComputeLogJointPsiMarginal() ;

    outputValue = pointedForest->GetLogJointPsiMarginal() ;
    // Rprintf("Marginal joint Psi: %.4e \n \n \n", pointedForest->GetLogJointPsiMarginal()) ;
  // }
  // else
  // {
  //   throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  // }
  return outputValue ;
}

// [[Rcpp::export]]

Eigen::VectorXd GetFullCondMean(SEXP forestPointer) {
  vec outputVec ;
  if (!(forestPointer == NULL))
  {
    XPtr<Forest> pointedForest(forestPointer) ; // Becomes a regular pointer again.
    outputVec = pointedForest->GetFullCondMean() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return outputVec ;
}

// [[Rcpp::export]]

Eigen::VectorXd GetFullCondSDs(SEXP forestPointer) {
  vec outputVec ;
  if (!(forestPointer == NULL))
  {
    XPtr<Forest> pointedForest(forestPointer) ; // Becomes a regular pointer again.
    outputVec = pointedForest->GetFullCondSDs() ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return outputVec ;
}

// [[Rcpp::export]]

Rcpp::List ComputeCondPredStats(SEXP forestPointer, NumericMatrix spCoordsForPredict, NumericVector timeForPredict,
                                 NumericMatrix covariateMatrixForPredict) {
  vec Evar, Hmean ;
  if (!(forestPointer == NULL))
  {
    XPtr<Forest> pointedForest(forestPointer) ; // Becomes a regular pointer again.

    // pointedForest->ComputeHpred() ;

    Hmean = pointedForest->GetHmatPred() * pointedForest->GetFullCondMean() ;
    Rcout << "Computing Evar..." << std::endl ;
    Evar = pointedForest->ComputeEvar() ;
    Rcout << "Done!" << std::endl ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null." ) ;
  }
  return Rcpp::List::create(Named("Hmean") = Hmean, Named("Evar") = Evar) ;
}

// [[Rcpp::export]]

int GetNumTips(SEXP forestPointer) {
  XPtr<AugTree> pointedForest(forestPointer) ; // Becomes a regular pointer again.
  int value = pointedForest->GetNumTips() ;
  return value ;
}

// [[Rcpp::export]]

Eigen::VectorXi GetPredObsOrder(SEXP forestPointer) {
  XPtr<Forest> pointedForest(forestPointer) ; // Becomes a regular pointer again.
  uvec value = pointedForest->GetObsOrderForFpredMat() ;
  return value ;
}

// [[Rcpp::export]]

Eigen::SparseMatrix<double> GetHmat(SEXP forestPointer) {
  XPtr<Forest> pointedForest(forestPointer) ; // Becomes a regular pointer again.
  Eigen::SparseMatrix<double> value = pointedForest->GetHmat() ;
  return value ;
}
