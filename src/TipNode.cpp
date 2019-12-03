#include "TipNode.h"

using namespace Eigen ;
using namespace MRAinla ;
using namespace std ;

void TipNode::computeUpred(const maternVec & covParasSp, const spatialcoor & predictLocations,
                           const double & spaceNuggetSD, const string & distMethod) {
  if (m_predsInNode.size() > 0) {
    std::vector<TreeNode *> brickList = getAncestors() ;

    spatialcoor subLocations = predictLocations.subset(m_predsInNode) ; // Will I need to invoke the destructor to avoid a memory leak?
    m_UmatList.at(0) = computeCovMat(subLocations, brickList.at(0)->GetKnotsCoor(), covParasSp, spaceNuggetSD, distMethod) ;

    for (uint l = 1; l <= m_depth; l++) {
      mat firstTerm = computeCovMat(subLocations, brickList.at(l)->GetKnotsCoor(), covParasSp, spaceNuggetSD, distMethod) ;
      mat secondTerm = mat::Zero(firstTerm.rows(), firstTerm.cols()) ;
      for (uint k = 0 ; k <= l-1; k++) {
         secondTerm += m_UmatList.at(k) * brickList.at(k)->GetKmatrix().selfadjointView<Upper>() * brickList.at(l)->GetWlist().at(k).transpose() ;
      }
      m_UmatList.at(l) = firstTerm - secondTerm ;
    }
  }
}

void TipNode::SetPredictLocations(const inputdata & predictData) {
  ArrayXi indices = deriveObsInNode(predictData) ;
  m_predsInNode = indices ;
}
