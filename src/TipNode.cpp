#include "TipNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void TipNode::DeriveOmega(const vec & responseValues) {
  vec subResponses = responseValues.elem(m_obsInNode) ;
  for (uint i = 0 ; i < m_depth; i++) {
    m_omegaTilde.at(i) = GetB(i).transpose() * m_SigmaInverse * subResponses ;
  }
}

void TipNode::computeUpred(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const spatialcoor & predictLocations,
                           const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {
  if (m_predsInNode.size() > 0) {
    std::vector<TreeNode *> brickList = getAncestors() ;
    m_UmatList.resize(m_depth + 1) ;

    spatialcoor subLocations = predictLocations.subset(m_predsInNode) ; // Will I need to invoke the destructor to avoid a memory leak?
    m_UmatList.at(0) = computeCovMat(subLocations, brickList.at(0)->GetKnotsCoor(), covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;

    for (uint l = 1; l <= m_depth; l++) {
      mat firstTerm = computeCovMat(subLocations, brickList.at(l)->GetKnotsCoor(), covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;
      mat secondTerm(firstTerm.n_rows, firstTerm.n_cols, fill::zeros) ;
      for (uint k = 0 ; k <= l-1; k++) {
         secondTerm += m_UmatList.at(k) * brickList.at(k)->GetKmatrix() * trans(brickList.at(l)->GetWlist().at(k)) ;
      }
      m_UmatList.at(l) = firstTerm - secondTerm ;
    }
  }
}

void TipNode::SetPredictLocations(const inputdata & predictData) {
  uvec indices = deriveObsInNode(predictData) ;
  m_predsInNode = indices ;
}
