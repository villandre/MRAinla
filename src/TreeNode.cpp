#include <math.h>

#include "TreeNode.h"

using namespace Eigen ;
using namespace MRAinla ;
using namespace std ;

ArrayXi TreeNode::deriveObsInNode(const spatiotempcoor & dataset) {
  Eigen::Array<bool, Dynamic, 1> lonCheck = (dataset.spatialCoords.col(0) > m_dimensions.longitude.matrix().minCoeff()) *
    (dataset.spatialCoords.col(0) <= m_dimensions.longitude.matrix().maxCoeff()) ; // Longitude check
  Eigen::Array<bool, Dynamic, 1> latCheck = (dataset.spatialCoords.col(1) > m_dimensions.latitude.matrix().minCoeff()) *
    (dataset.spatialCoords.col(1) <= m_dimensions.latitude.matrix().maxCoeff()) ; // Latitude check
  Eigen::Array<bool, Dynamic, 1> timeCheck(dataset.timeCoords)
  double timeValue = m_dimensions.timeValue ;
  for (uint i = 0; i < dataset.timeCoords.size(); i++) {
    timeCheck(i) = (dataset.timeCoords(i) == timeValue) ;
  }
  ArrayXi output = find(lonCheck * latCheck * timeCheck) ;
  return output ;
}

double TreeNode::MaternCovFunction(const double & distance, const maternVec & covParasSpace,
                                   const double & spaceNuggetSD) {
  double spExp = maternCov(distance, covParasSpace.m_rho, covParasSpace.m_smoothness , covParasSpace.m_scale,
                           spaceNuggetSD) ;

  return spExp  ;
}

void TreeNode::baseComputeWmat(const maternVec & covParasSp, const double & spaceNuggetSD, const string & distMethod) {

  std::vector<TreeNode *> brickList = getAncestors() ;

  m_Wlist.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor(), covParasSp, spaceNuggetSD, distMethod) ;

  for (uint l = 1; l <= m_depth; l++) {
    mat firstMat = computeCovMat(m_knotsCoor, brickList.at(l)->GetKnotsCoor(), covParasSp, spaceNuggetSD, distMethod) ;
    mat secondMat = mat::Zero(firstMat.rows(), firstMat.cols()) ;
    for (uint k = 0; k < l ; k++) {
      // if (k < m_depth) {
      mat increment = (m_Wlist.at(k) *
        brickList.at(k)->GetKmatrix().selfadjointView<Upper>() *
        brickList.at(l)->GetWlist().at(k).transpose()) ;
      secondMat += increment ;
      // } else {
      //   secondMat.noalias() += m_Wlist.at(m_depth) *
      //     brickList.at(k)->GetKmatrix().selfadjointView<Upper>() *
      //     brickList.at(l)->GetWlist().at(m_depth).transpose() ;
      // }
    }

    m_Wlist.at(l) = firstMat - secondMat ;
  }
}

std::vector<TreeNode *> TreeNode::getAncestors() {
  TreeNode * currentAddress = this ;
  std::vector<TreeNode *> ancestorsList ;
  ancestorsList.push_back(currentAddress) ;
  if (m_parent == currentAddress) {
    return ancestorsList ;
  }
  for(uint i = 0; i < m_depth; i++) {
    currentAddress = currentAddress->GetParent() ;
    ancestorsList.push_back(currentAddress) ;
  }
  std::reverse(ancestorsList.begin(), ancestorsList.end()) ;
  return ancestorsList ;
}

mat TreeNode::computeCovMat(const spatialcoor & sp1, const spatialcoor & sp2, const maternVec & covParasSp, const double & spaceNuggetSD, const string & distMethod) {

  mat covMat = mat::Zero(sp1.spatialCoords.rows(), sp2.spatialCoords.rows()) ;
  for (uint rowIndex = 0; rowIndex < sp1.spatialCoords.rows() ; rowIndex++) {
    for (uint colIndex = 0; colIndex < sp2.spatialCoords.rows() ; colIndex++) {
      ArrayXd space1 = sp1.spatialCoords.row(rowIndex) ;
      ArrayXd space2 = sp2.spatialCoords.row(colIndex) ;

      double rangeValue = spDistance(space1, space2, distMethod) ;
      covMat(rowIndex, colIndex) = MaternCovFunction(rangeValue, covParasSp, spaceNuggetSD) ;
    }
  }
  return covMat ;
}
