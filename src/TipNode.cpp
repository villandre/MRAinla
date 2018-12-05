#include "TipNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void TipNode::ComputeParasEtaDeltaTilde(const spatialcoor & predictLocations, const inputdata & dataset, const arma::vec & fixedParas, const arma::vec & covParas) {
  computeUpred(covParas, predictLocations) ;
  computeVpred(covParas, predictLocations) ;

  computeDeltaTildeParas(dataset) ;
}

void TipNode::computeUpred(const vec & covParas, const spatialcoor & predictLocations) {
  if (m_predictLocIndices.size() > 0) {
    std::vector<TreeNode *> brickList = getAncestors() ;
    m_UmatList.resize(m_depth + 1) ;

    spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ; // Will I need to invoke the destructor to avoid a memory leak?
    m_UmatList.at(0) = computeCovMat(subLocations, m_knotsCoor, covParas) ;

    mat firstTerm(m_UmatList.at(0).n_rows, m_UmatList.at(0).n_cols) ;
    mat secondTerm(firstTerm.n_rows, firstTerm.n_cols) ;
    for (uint l = 1; l <= m_depth; l++) {
      firstTerm = computeCovMat(subLocations, brickList.at(l)->GetKnotsCoor(), covParas) ;
      secondTerm.fill(0) ;
      for (uint k = 0 ; k <= l-1; k++) {
         secondTerm += m_UmatList.at(k) * brickList.at(k)->GetKmatrix() * trans(brickList.at(l)->GetWlist().at(k)) ;
      }
      m_UmatList.at(l) = firstTerm - secondTerm ;
    }
  }
}

void TipNode::computeVpred(const arma::vec & covParas, const spatialcoor & predictLocations) {
  if (m_predictLocIndices.size() > 0) {
    std::vector<TreeNode *> brickList = getAncestors() ;
    spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ;
    mat firstTerm = computeCovMat(subLocations, subLocations, covParas) ;
    mat secondTerm(firstTerm.n_rows, firstTerm.n_cols) ;
    for (uint k = 0 ; k <= m_depth-1; k++) {
      secondTerm += m_UmatList.at(k) * brickList.at(k)->GetKmatrix() * trans(m_UmatList.at(k)) ;
    }
    m_V = firstTerm - secondTerm ;
  }
};

void TipNode::computeDeltaTildeParas(const inputdata & dataset, const vec & fixedEffValues) {
  arma::vec subResponses = dataset.responseValues.elem(m_obsInNode) ;
  arma::vec fixedSubtraction(m_obsInNode.size(), 0) ;
  arma::mat incCovariates = dataset.covariateValues ;
  incCovariates.insert_cols(0, 1) ;
  incCovariates.col(0).fill(1) ;
  fixedSubtraction = incCovariates.rows(m_obsInNode) * fixedEffValues ;
  subResponses -= fixedSubtraction ;
  m_deltaTilde.meanPara = m_UmatList.at(m_depth) * m_SigmaInverse * subResponses ;
  m_deltaTilde.covPara = m_V - GetLM() * m_SigmaInverse * trans(GetLM()) ;
}

void TipNode::deriveBtilde(const spatialcoor & predictLocations) {
  m_Btilde.resize(m_depth + 1) ;
  for (auto & i : m_Btilde) {
    i.resize(m_depth) ;
  }
  mat bPredMat = computeBpred(predictLocations) ;
  m_Btilde.at(m_depth).at(m_depth - 1) = bPredMat - GetLM() * m_SigmaInverse * GetB(m_depth - 1) ;

};

// .computeBpreds <- function(gridObj, locations) {
//
// allTips <- .tipAddresses(gridObj)
//
//   recurseBpred <- function(tipAddress) {
//     subLocations <- subset(x = locations, latExtent = tipAddress$dimensions$latitude, lonExtent = tipAddress$dimensions$longitude, timeExtent = tipAddress$dimensions$time)
//     tipAddress$bPred <- vector('list', length = gridObj$M+1)
//     if (length(subLocations) == 0) {
//       tipAddress$bPred <- NULL
//       return(invisible())
//     }
//     tipAddress$bPred[[1]] <- gridObj$covFct(subLocations, gridObj$knotPositions)
//
//       brickList <- .getAllParentAddresses(brickObj = tipAddress, flip = FALSE) ## Lists bricks from root to current tip (included).
//
//     recurseBpred <- function(level) {
//       currentV <- gridObj$covFct(subLocations, brickList[[level+1]]$knotPositions) - tipAddress$bPred[[1]] %*% gridObj$K %*% t(brickList[[level+1]]$bKnots[[1]]) ### SAME ISSUE AS IN .computeBknots with the position of the t(...). ###
//       if (level == 1) {
//         tipAddress$bPred[[2]] <- currentV
//         return(invisible())
//       }
//
//       for (i in 2:level) {
//         currentV <- currentV - tipAddress$bPred[[i]] %*% brickList[[i]]$K %*% t(brickList[[level+1]]$bKnots[[i]]) # Changed the position of the t to ensure correspondence.
//       }
//       tipAddress$bPred[[level+1]] <- currentV
//         invisible()
//     }
//
//     for (levelIndex in 1:gridObj$M) {
//       recurseBpred(levelIndex)
//     }
//     invisible()
//   }
//   lapply(allTips, recurseBpred)
// }

