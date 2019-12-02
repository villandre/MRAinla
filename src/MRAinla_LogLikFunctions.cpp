#include <iostream>
#include <string>
#include <algorithm>
#include <RcppEigen.h>
// #include "gperftools/profiler.h"

#include "AugTree.h"

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
  std::vector<double> timeAsVector = Rcpp::as<std::vector<double>>(obsTime) ;
  std::sort(timeAsVector.begin(), timeAsVector.end()) ;
  timeAsVector.erase(unique( timeAsVector.begin(), timeAsVector.end() ),
                      timeAsVector.end() );
  std::vector<AugTree *> MRAgridVector ;
  for (auto & i : timeAsVector) {
    Eigen::Array<bool, Eigen::Dynamic, 1> logicalTest = (time == i) ;
    Eigen::ArrayXi timeIndices = find(logicalTest) ;
    Eigen::ArrayXXd subObsSp = rows(sp, timeIndices) ;
    Eigen::ArrayXXd subCovar = rows(covariateMat, timeIndices) ;
    vec subResponse = elem(responseArray, timeIndices).matrix() ;
    ArrayXd subTime = elem(time, timeIndices) ;

    AugTree * MRAgrid = new AugTree(Mlon, Mlat, lonR, latR, subResponse, subObsSp, subTime, predCovariates, predSp, predTimeVec, seedForRNG, subCovar, numKnotsRes0, J, dMethod) ;
    MRAgridVector.push_back(MRAgrid) ;
  }

  XPtr<std::vector<AugTree *>> p(&MRAgridVector, false) ; // Disabled automatic garbage collection.

  return List::create(Named("gridPointer") = p) ;
}

// [[Rcpp::export]]

double LogJointHyperMarginalToWrap(SEXP treePointer, Rcpp::List MRAhyperparas,
         double timeCovPara, double fixedEffSD, double errorSD, Rcpp::List MRAcovParasGammaAlphaBeta,
         Rcpp::NumericVector FEmuVec, NumericVector fixedEffGammaAlphaBeta,
         NumericVector errorGammaAlphaBeta, NumericVector timeGammaAlphaBeta,
         double spaceNuggetSD, bool recordFullConditional, bool processPredictions) {
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
      pointedTree->SetTimeCovParaGammaAlphaBeta(GammaHyperParas(timeGammaAlphaBeta(0), timeGammaAlphaBeta(1))) ;
      pointedTree->SetFEmu(FEmu) ;
      pointedTree->SetSpaceNuggetSD(spaceNuggetSD) ;
      std::vector<TreeNode *> tipNodes = pointedTree->GetLevelNodes(pointedTree->GetM()) ;
      for (auto & i : tipNodes) {
        i->SetUncorrSD(0.001) ; // Is this a nugget effect?
      }
    }
    pointedTree->SetErrorSD(errorSD) ;
    pointedTree->SetFixedEffSD(fixedEffSD) ;

    pointedTree->SetMRAcovParas(MRAhyperparas) ;
    pointedTree->SetRecordFullConditional(recordFullConditional) ;
    pointedTree->SetProcessPredictions(processPredictions) ;

    pointedTree->ComputeLogJointPsiMarginal() ;

    outputValue = pointedTree->GetLogJointPsiMarginal() ;
    // Rprintf("Marginal joint Psi: %.4e \n \n \n", pointedTree->GetLogJointPsiMarginal()) ;
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
                                 NumericMatrix covariateMatrixForPredict) {
  vec Evar, Hmean ;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.

    // pointedTree->ComputeHpred() ;

    Hmean = pointedTree->GetHmatPred() * pointedTree->GetFullCondMean() ;
    Rcout << "Computing Evar..." << std::endl ;
    Evar = pointedTree->ComputeEvar() ;
    Rcout << "Done!" << std::endl ;
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

// [[Rcpp::export]]

Eigen::SparseMatrix<double> GetHmat(SEXP treePointer) {
  XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.
  Eigen::SparseMatrix<double> value = pointedTree->GetHmat() ;
  return value ;
}
