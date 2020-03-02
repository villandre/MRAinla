#include <iostream>
#include <string>
#include <algorithm>
#include <RcppEigen.h>
// #include <unsupported/Eigen/SparseExtra>
// #include "gperftools/profiler.h"

#include "AugTree.h"

using namespace Eigen;
using namespace Rcpp;
using namespace MRAinla;
using namespace std;

// [[Rcpp::export]]

List setupNestedGrids(NumericVector responseValues,
                      NumericMatrix spCoords,
                      NumericMatrix predCoords,
                      NumericVector obsTime,
                      NumericVector predTime,
                      NumericMatrix covariateMatrix,
                      NumericMatrix predCovariateMatrix,
                      uint Mlon, uint Mlat, uint Mtime,
                      NumericVector lonRange, NumericVector latRange, NumericVector timeRange,
                      uint randomSeed,
                      int numKnotsRes0,
                      double J,
                      String distMethod,
                      Rcpp::List MaternParsHyperpars,
                      Rcpp::NumericVector fixedEffParsHyperpars,
                      NumericVector errorParsHyperpars,
                      Rcpp::NumericVector FEmuVec,
                      double nuggetSD,
                      bool normalHyperprior,
                      double tipKnotsThinningRate)
{
  ArrayXd lonRinit = as<ArrayXd>(lonRange) ;
  ArrayXd latRinit = as<ArrayXd>(latRange) ;
  ArrayXd timeRinit = as<ArrayXd>(timeRange) ;
  Array2d lonR = lonRinit.segment(0,2) ;
  Array2d latR = latRinit.segment(0,2) ;
  Array2d timeR = timeRinit.segment(0,2) ;
  vec response = as<vec>(responseValues) ;
  ArrayXXd sp = as<ArrayXXd>(spCoords) ;
  ArrayXXd predSp, predCovariates ;
  ArrayXd predTimeVec ;
  if (!(predCoords == R_NilValue)) {
    ArrayXXd predSp = as<ArrayXXd>(predCoords) ;
    ArrayXd predTimeVec = as<ArrayXd>(predTime) ;
    ArrayXXd predCovariates = as<ArrayXXd>(predCovariateMatrix) ;
  }
  ArrayXd time = as<ArrayXd>(obsTime) ;
  string dMethod = as<std::string>(Rcpp::wrap(distMethod)) ;

  unsigned long int seedForRNG = randomSeed ;

  ArrayXXd covariateMat = as<ArrayXXd>(covariateMatrix) ;

  AugTree * MRAgrid = new AugTree(Mlon, Mlat, Mtime,
                                  lonR, latR, timeR,
                                  response, sp, time,
                                  predCovariates, predSp, predTimeVec,
                                  seedForRNG,
                                  covariateMat,
                                  numKnotsRes0, J,
                                  dMethod,
                                  MaternParsHyperpars, fixedEffParsHyperpars, errorParsHyperpars,
                                  FEmuVec,
                                  nuggetSD,
                                  normalHyperprior,
                                  tipKnotsThinningRate) ;
  XPtr<AugTree> p(MRAgrid, false) ; // Disabled automatic garbage collection.
  return List::create(Named("nestedGridsPointer") = p) ;
}

// [[Rcpp::export]]

double LogJointHyperMarginalToWrap(SEXP treePointer, Rcpp::List MaternHyperpars,
         double fixedEffSD, double errorSD, bool recordFullConditional,
         bool processPredictions) {
  mat posteriorMatrix ;
  double outputValue = 0 ;

  XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.

  std::vector<TreeNode *> tipNodes = pointedTree->GetLevelNodes(pointedTree->GetM()) ;

  for (auto & i : tipNodes) {
    i->SetUncorrSD(0.001) ; // Is this a nugget effect?
  }
  pointedTree->SetMaternPars(MaternHyperpars) ;
  pointedTree->SetErrorSD(errorSD) ;
  pointedTree->SetFixedEffSD(fixedEffSD) ;

  pointedTree->SetRecordFullConditional(recordFullConditional) ;
  pointedTree->SetProcessPredictions(processPredictions) ;

  pointedTree->ComputeLogJointPsiMarginal() ;

  outputValue = pointedTree->GetLogJointPsiMarginal() ;
  // Rprintf("Marginal joint Psi: %.4e \n \n \n", pointedTree->GetLogJointPsiMarginal()) ;

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

Rcpp::List ComputeCondPredStats(SEXP treePointer) {
  vec Evar, Hmean ;
  if (!(treePointer == NULL))
  {
    XPtr<AugTree> pointedTree(treePointer) ; // Becomes a regular pointer again.

    Hmean = pointedTree->GetHmatPred() * pointedTree->GetFullCondMean() ;
    Rcout << "Computing Evar..." << std::endl ;
    Evar = pointedTree->ComputeEvar() ;
    Rcout << "Done!" << std::endl ;
  }
  else
  {
    throw Rcpp::exception("Pointer to MRA grid is null.") ;
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
