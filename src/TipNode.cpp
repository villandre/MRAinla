#include "TipNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void TipNode::DeriveOmega(const arma::vec & responseValues) {
  arma::vec subResponses = responseValues.elem(m_obsInNode) ;
  for (uint i = 0 ; i < m_depth; i++) {
    m_omegaTilde.at(i) = arma::trans(GetB(i)) * m_SigmaInverse * subResponses ;
  }
}

// void TipNode::ComputeParasEtaDeltaTilde(const spatialcoor & predictLocations, const inputdata & dataset, const arma::vec & covParas) {
//   if (m_predictLocIndices.size() > 0) {
//     computeUpred(covParas, predictLocations) ;
//     computeVpred(covParas, predictLocations) ;
//     computeDeltaTildeParas(dataset) ;
//   }
// }

// void TipNode::computeUpred(const vec & covParas, const spatialcoor & predictLocations) {
//   if (m_predictLocIndices.size() > 0) {
//     std::vector<TreeNode *> brickList = getAncestors() ;
//     m_UmatList.resize(m_depth + 1) ;
//
//     spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ; // Will I need to invoke the destructor to avoid a memory leak?
//     m_UmatList.at(0) = computeCovMat(subLocations, brickList.at(0)->GetKnotsCoor(), covParas) ;
//
//     for (uint l = 1; l <= m_depth; l++) {
//       mat firstTerm = computeCovMat(subLocations, brickList.at(l)->GetKnotsCoor(), covParas) ;
//       mat secondTerm(firstTerm.n_rows, firstTerm.n_cols, fill::zeros) ;
//       for (uint k = 0 ; k <= l-1; k++) {
//          secondTerm += m_UmatList.at(k) * brickList.at(k)->GetKmatrix() * trans(brickList.at(l)->GetWlist().at(k)) ;
//       }
//       m_UmatList.at(l) = firstTerm - secondTerm ;
//     }
//   }
// }

// void TipNode::computeVpred(const arma::vec & covParas, const spatialcoor & predictLocations) {
//   if (m_predictLocIndices.size() > 0) {
//     std::vector<TreeNode *> brickList = getAncestors() ;
//     spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ;
//     mat firstTerm = computeCovMat(subLocations, subLocations, covParas) ;
//     mat secondTerm = zeros<mat>(firstTerm.n_rows, firstTerm.n_cols) ;
//     for (uint k = 0 ; k <= m_depth-1; k++) {
//       secondTerm += m_UmatList.at(k) * brickList.at(k)->GetKmatrix() * trans(m_UmatList.at(k)) ;
//     }
//     m_V = firstTerm - secondTerm ;
//   }
// };

// void TipNode::computeDeltaTildeParas(const inputdata & dataset) {
//   if (m_predictLocIndices.size() > 0) {
//     arma::vec subResponses = dataset.responseValues.elem(m_obsInNode) ;
//     m_deltaTilde.meanPara = m_UmatList.at(m_depth) * m_SigmaInverse * subResponses ;
//     m_deltaTilde.covPara = m_V - GetLM() * m_SigmaInverse * trans(GetLM()) ;
//   }
// }
//
// void TipNode::deriveBtilde(const spatialcoor & predictLocations) {
//   if (m_predictLocIndices.size() > 0) {
//     spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ;
//
//     m_Btilde.resize(m_depth + 1) ;
//     for (uint i = 1; i < (m_depth + 1) ; i++) {
//       m_Btilde.at(i).resize(i) ;
//     }
//     std::vector<TreeNode *> brickList = getAncestors() ;
//     for (uint k = 0 ; k < m_depth ; k++) {
//       m_Btilde.at(m_depth).at(k) = m_bPred.at(k) - GetLM() * m_SigmaInverse * m_Wlist.at(k) ;
//     }
//
//     for (uint k = m_depth-1 ; k >= 1 ; k--) {
//       recurseBtilde(k , k-1) ;
//       if (k > 1) recurseBtilde(k, k-2) ;
//       if (k > 2) recurseBtilde(k, k -3) ;
//     }
//   }
// };
//
// void TipNode::recurseBtilde(const uint l, const uint k) {
//   std::vector<TreeNode *> brickList = getAncestors() ;
//   m_Btilde.at(l).at(k) = m_Btilde.at(l+1).at(k) - m_Btilde.at(l+1).at(l) * brickList.at(l)->GetKtilde() * brickList.at(l)->GetAlist().at(l).at(k) ;
// }

// void TipNode::computeBpred(const spatialcoor & predictLocations, const vec & covParas) {
//   if (m_predictLocIndices.size() > 0) {
//     spatialcoor subLocations = predictLocations.subset(m_predictLocIndices) ;
//     m_bPred.resize(m_depth+1) ;
//     std::vector<TreeNode *> brickList = getAncestors() ;
//     m_bPred.at(0) = computeCovMat(subLocations, brickList.at(0)->GetKnotsCoor(), covParas) ;
//
//     for (uint levelIndex = 1 ; levelIndex <= m_depth ; levelIndex++) {
//       mat currentV = computeCovMat(subLocations, brickList.at(levelIndex)->GetKnotsCoor(), covParas) -
//         m_bPred.at(0) * brickList.at(0)->GetKmatrix() * trans(brickList.at(levelIndex)->GetBknots().at(0)) ;
//       if (levelIndex > 1) {
//         for (uint i = 1 ; i < levelIndex ; i++) {
//           currentV -= m_bPred.at(i) * brickList.at(i)->GetKmatrix() * trans(brickList.at(levelIndex)->GetBknots().at(i)) ;
//         }
//       }
//       m_bPred.at(levelIndex) = currentV ;
//     }
//   }
// }

// GaussDistParas TipNode::CombineEtaDelta(const inputdata & dataset, const vec & fixedEffParas) {
//   vec defaultVec = zeros<vec>(1) ;
//   mat defaultMat = zeros<mat>(1,1) ;
//   defaultMat.at(0,0) = -1 ;
//   GaussDistParas estimatesFromRegion(defaultVec, defaultMat) ;
//
//   if (m_predictLocIndices.size() > 0) {
//
//     vec meanVec = zeros<vec>(m_predictLocIndices.size()) ;
//     mat covMat = zeros<mat>(m_predictLocIndices.size(), m_predictLocIndices.size()) ;
//     estimatesFromRegion.meanPara.resize(meanVec.size()) ;
//     estimatesFromRegion.covPara.resize(covMat.n_rows, covMat.n_cols) ;
//     std::vector<TreeNode *> brickList = getAncestors() ;
//     for (uint m=0 ; m <= m_depth-1 ; m++) {
//       meanVec += m_Btilde.at(m+1).at(m) * brickList.at(m)->GetEtaDelta().meanPara ;
//       covMat += m_Btilde.at(m+1).at(m) * brickList.at(m)->GetEtaDelta().covPara * trans(m_Btilde.at(m+1).at(m)) ;
//     }
//
//     meanVec += m_deltaTilde.meanPara ;
//
//     // We need to re-add the mean contribution of the fixed effect parameters, which were subtracted in the beginning with AugTree::CenterResponse
//
//     mat subCovar = dataset.covariateValues.rows(m_predictLocIndices) ;
//     vec interceptVec = ones<vec>(m_predictLocIndices.size()) ;
//     subCovar.insert_cols(0, interceptVec) ;
//     meanVec += subCovar * fixedEffParas ;
//     //
//     covMat += m_deltaTilde.covPara ;
//     estimatesFromRegion.meanPara = meanVec ;
//     estimatesFromRegion.covPara = covMat ;
//   }
//   return estimatesFromRegion ;
// }
