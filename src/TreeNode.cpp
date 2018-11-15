#include <gsl/gsl_sf_exp.h>

#include "TreeNode.h"

using namespace arma ;
using namespace MRAinla ;

void TreeNode::deriveObsInNode(const inputdata & dataset) {
  uvec lonCheck = (dataset.spatialCoords.col(0) > min(m_dimensions.longitude)) %
    (dataset.spatialCoords.col(0) <= max(m_dimensions.longitude)) ; // Longitude check
  uvec latCheck = (dataset.spatialCoords.col(1) > min(m_dimensions.latitude)) %
    (dataset.spatialCoords.col(1) <= max(m_dimensions.latitude)) ; // Latitude check
  uvec timeCheck = (dataset.timeCoords > min(m_dimensions.time)) %
    (dataset.timeCoords <= max(m_dimensions.time)) ; // Time check
  m_obsInNode = find(lonCheck % latCheck % timeCheck) ; // find is equivalent to which in R
}

double TreeNode::covFunction(const Spatiotemprange & distance) {
  double spExp = -(distance.sp/m_covPara.at(0)) ;
  double timeExp = -(distance.time/m_covPara.at(1)) ;
  return gsl_sf_exp(spExp) * gsl_sf_exp(timeExp) ;
};

mat TreeNode::ComputeBaseKmat() {
  mat covMatrix(m_knotsCoor.timeCoords.size(), m_knotsCoor.timeCoords.size(), fill::zeros) ;
  Spatiotemprange container ;
  for (uint rowIndex = 0 ; rowIndex < m_knotsCoor.timeCoords.size() ; rowIndex++) {
    for (uint colIndex = 0 ; colIndex <= rowIndex ; colIndex++) {
       container = sptimeDistance(m_knotsCoor.spatialCoords.row(rowIndex),
                                  m_knotsCoor.timeCoords.at(rowIndex),
                                  m_knotsCoor.spatialCoords.row(colIndex),
                                  m_knotsCoor.timeCoords.at(colIndex)) ;
      covMatrix.at(rowIndex, colIndex) = covFunction(container) ;
    }
  }
  return symmatl(covMatrix) ;
}

void TreeNode::ComputeWmat() {
  uint M = GetM() ;
  if (m_depth == 0) {
    m_Wlist.at(0) = m_Kinverse ;
  } else {
    std::vector<TreeNode *> brickList = getAncestors() ;
    m_Wlist.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor()) ;

    for (uint l = 1; l < m_depth+1; l++) {
      mat firstMat = computeCovMat(m_knotsCoor, brickList.at(l)->GetKnotsCoor()) ;
      mat secondMat(firstMat.n_rows, firstMat.n_cols, fill::zeros) ;
      for (uint k = 0; k < l ; k++) {
        secondMat = secondMat + m_Wlist.at(k) *
          brickList.at(k)->GetKmatrix() *
          trans(brickList.at(l)->GetWlist().at(k)) ;
      }
      m_Wlist.at(l) = firstMat - secondMat ;
    }
    if (m_depth < M) {
      m_Kinverse = m_Wlist.at(m_depth+1) ;
      m_K = inv_sympd(m_Kinverse) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
    }
  }
}

std::vector<TreeNode *> TreeNode::getAncestors() {
  TreeNode * currentAddress = this ;
  std::vector<TreeNode *> ancestorsList(m_depth+1) ;
  ancestorsList.at(m_depth) = currentAddress ;
  if (m_parent == currentAddress) {
    return ancestorsList ;
  }
  for(uint i = 0; i < m_depth; i++) {
    currentAddress = currentAddress->GetParent() ;
    ancestorsList.at(currentAddress->GetDepth()) = currentAddress ;
  }
  return ancestorsList ;
}

mat TreeNode::computeCovMat(const spatialcoor & spTime1, const spatialcoor & spTime2) {

  mat covMat(spTime1.timeCoords.size(), spTime2.timeCoords.size(), fill::zeros) ;
  for (uint rowIndex = 0; rowIndex < spTime1.timeCoords.size() ; rowIndex++) {
    for (uint colIndex = 0; colIndex < spTime2.timeCoords.size() ; colIndex++) {
      Spatiotemprange rangeValue = sptimeDistance(spTime1.spatialCoords.row(rowIndex),
                                                  spTime1.timeCoords.at(rowIndex),
                                                  spTime2.spatialCoords.row(colIndex),
                                                  spTime2.timeCoords.at(colIndex)) ;
      covMat.at(rowIndex, colIndex) = covFunction(rangeValue) ;
    }
  }
  return symmatl(covMat) ;
}
