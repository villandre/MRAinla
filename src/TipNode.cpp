#include "TipNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void TipNode::ComputeParasEtaDeltaTilde(const spatialcoor & predictLocations, const inputdata & dataset, const arma::vec & fixedParas, const arma::vec & covParas) {
  computeUpred(covParas, predictLocations) ;
  computeVpred(covParas, predictLocations) ;

  computeDeltaTildeParas(dataset, fixedParas) ;
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
  if (m_predictLocIndices.size() > 0) {
    spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ;

    m_Btilde.resize(m_depth + 1) ;
    for (uint i = 1; i < (m_depth + 1) ; i++) {
      m_Btilde.at(i).resize(i) ;
    }
    std::vector<TreeNode *> brickList = getAncestors() ;
    for (uint k = 0 ; k < m_depth ; k++) {
      m_Btilde.at(m_depth).at(k) = m_bPred.at(k) - GetLM() * m_SigmaInverse * m_Wlist.at(k) ;
    }

    for (uint k = m_depth-1 ; k >= 1 ; k--) {
      recurseBtilde(k , k-1) ;
      if (k > 1) recurseBtilde(k, k-2) ;
      if (k > 2) recurseBtilde(k, k -3) ;
    }
  }
};

void TipNode::recurseBtilde(const uint l, const uint k) {
  std::vector<TreeNode *> brickList = getAncestors() ;
  m_Btilde.at(l).at(k) = m_Btilde.at(l+1).at(k) - m_Btilde.at(l+1).at(l) * brickList.at(l)->GetKtilde() * brickList.at(l)->GetAlist().at(l).at(k) ;
}

void TipNode::computeBpred(const spatialcoor & predictLocations, const vec & covParas) {
  if (m_predictLocIndices.size() > 0) {
    spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ;
    m_bPred.resize(m_depth+1) ;
    m_bPred.at(0) = computeCovMat(subLocations, m_knotsCoor, covParas) ;
    std::vector<TreeNode *> brickList = getAncestors() ;

    for (uint levelIndex = 1 ; levelIndex <= m_depth ; levelIndex++) {
      mat currentV = computeCovMat(subLocations, brickList.at(levelIndex)->GetKnotsCoor(), covParas) -
        m_bPred.at(0) * brickList.at(0)->GetKmatrix() * trans(brickList.at(levelIndex)->GetBknots().at(0)) ;
      if (levelIndex > 1) {
        for (uint i = 1 ; i < levelIndex ; i++) {
          currentV -= m_bPred.at(i) * brickList.at(i)->GetKmatrix() * trans(brickList.at(levelIndex)->GetBknots().at(i)) ;
        }
      }
      m_bPred.at(levelIndex) = currentV ;
    }
  }
}

GaussDistParas TipNode::CombineEtaDelta() {
  vec defaultVec(1,0) ;
  mat defaultMat(1,1,fill::zeros) ;
  GaussDistParas estimatesFromRegion(defaultVec, defaultMat) ;

  if (m_predictLocIndices.size() > 0) {

    vec meanVec(m_predictLocIndices.size(), 0) ;
    mat covMat(m_predictLocIndices.size(), m_predictLocIndices.size(), fill::zeros) ;
    std::vector<TreeNode *> brickList = getAncestors() ;
    for (uint m=0 ; m <= m_depth-1 ; m++) {
      meanVec += m_Btilde.at(m+1).at(m) * brickList.at(m)->GetEtaDelta().meanPara ;
      covMat += m_Btilde.at(m+1).at(m) * brickList.at(m)->GetEtaDelta().covPara * trans(m_Btilde.at(m+1).at(m)) ;
    }
    meanVec += m_deltaTilde.meanPara ;
    covMat += m_deltaTilde.covPara ;
    estimatesFromRegion.meanPara = meanVec ;
    estimatesFromRegion.covPara = covMat ;
  }
  return estimatesFromRegion ;
}
