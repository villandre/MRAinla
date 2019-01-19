#include <gsl/gsl_sf_exp.h>

#include "TreeNode.h"

using namespace arma ;
using namespace MRAinla ;

arma::uvec TreeNode::deriveObsInNode(const spatialcoor & dataset) {
  uvec lonCheck = (dataset.spatialCoords.col(0) > min(m_dimensions.longitude)) %
    (dataset.spatialCoords.col(0) <= max(m_dimensions.longitude)) ; // Longitude check
  uvec latCheck = (dataset.spatialCoords.col(1) > min(m_dimensions.latitude)) %
    (dataset.spatialCoords.col(1) <= max(m_dimensions.latitude)) ; // Latitude check
  uvec timeCheck = (dataset.timeCoords > min(m_dimensions.time)) %
    (dataset.timeCoords <= max(m_dimensions.time)) ; // Time check
  return find(lonCheck % latCheck % timeCheck) ; // find is equivalent to which in R
}

double TreeNode::covFunction(const Spatiotemprange & distance, const vec & covParameters) {
  double spExp = -(distance.sp/covParameters.at(0)) ;
  double timeExp = -(distance.time/covParameters.at(1)) ;
  return gsl_sf_exp(spExp) * gsl_sf_exp(timeExp) ;
};

arma::mat TreeNode::ComputeCovMatrix(const vec & covParaVec) {
  Spatiotemprange container ;
  mat covMatrix(m_knotsCoor.timeCoords.size(), m_knotsCoor.timeCoords.size(), fill::zeros) ;
  for (uint rowIndex = 0 ; rowIndex < m_knotsCoor.timeCoords.size() ; rowIndex++) {
    for (uint colIndex = 0 ; colIndex <= rowIndex ; colIndex++) {
      vec space1 = conv_to<vec>::from(m_knotsCoor.spatialCoords.row(rowIndex)) ;
      uint time1 = m_knotsCoor.timeCoords.at(rowIndex) ;
      vec space2 = conv_to<vec>::from(m_knotsCoor.spatialCoords.row(colIndex)) ;
      uint time2 = m_knotsCoor.timeCoords.at(colIndex) ;
      container = sptimeDistance(space1, time1, space2, time2) ;
      covMatrix.at(rowIndex, colIndex) = covFunction(container, covParaVec) ;
    }
  }
  return covMatrix ;
}

void TreeNode::baseComputeWmat(const vec & covParas) {
  uint M = GetM() ;
  std::vector<TreeNode *> brickList = getAncestors() ;
  m_Wlist.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor(), covParas) ;

  for (uint l = 1; l < m_depth+1; l++) {
    mat firstMat = computeCovMat(m_knotsCoor, brickList.at(l)->GetKnotsCoor(), covParas) ;
    mat secondMat(firstMat.n_rows, firstMat.n_cols, fill::zeros) ;
    for (uint k = 0; k < l ; k++) {
      secondMat += m_Wlist.at(k) *
        brickList.at(k)->GetKmatrix() *
        trans(brickList.at(l)->GetWlist().at(k)) ;
    }
    m_Wlist.at(l) = firstMat - secondMat ;
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

mat TreeNode::computeCovMat(const spatialcoor & spTime1, const spatialcoor & spTime2, const vec & covParas) {
  mat covMat(spTime1.timeCoords.size(), spTime2.timeCoords.size(), fill::zeros) ;
  for (uint rowIndex = 0; rowIndex < spTime1.timeCoords.size() ; rowIndex++) {
    for (uint colIndex = 0; colIndex < spTime2.timeCoords.size() ; colIndex++) {
      vec space1 = conv_to<vec>::from(spTime1.spatialCoords.row(rowIndex)) ;
      uint time1 = spTime1.timeCoords.at(rowIndex) ;
      vec space2 = conv_to<vec>::from(spTime2.spatialCoords.row(colIndex)) ;
      uint time2 = spTime2.timeCoords.at(colIndex) ;
      Spatiotemprange rangeValue = sptimeDistance(space1, time1, space2, time2) ;
      covMat.at(rowIndex, colIndex) = covFunction(rangeValue, covParas) ;
    }
  }

  return covMat ;
}

void TreeNode::SetPredictLocations(const spatialcoor & predictLocations) {
  uvec lonCheck = (predictLocations.spatialCoords.col(0) > min(m_dimensions.longitude)) %
    (predictLocations.spatialCoords.col(0) <= max(m_dimensions.longitude)) ; // Longitude check
  uvec latCheck = (predictLocations.spatialCoords.col(1) > min(m_dimensions.latitude)) %
    (predictLocations.spatialCoords.col(1) <= max(m_dimensions.latitude)) ; // Latitude check
  uvec timeCheck = (predictLocations.timeCoords > min(m_dimensions.time)) %
    (predictLocations.timeCoords <= max(m_dimensions.time)) ; // Time check
  m_predictLocIndices = find(lonCheck % latCheck % timeCheck) ; // find is equivalent to which in R
}

void TreeNode::initiateBknots(const vec & covParas) {
  std::vector<TreeNode *> brickList = getAncestors() ;
  if (m_depth > 0) {
    m_bKnots.resize(m_depth) ;
    m_bKnots.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor(), covParas) ;
  }
}

void TreeNode::completeBknots(const vec & covParas, const uint level) {
  std::vector<TreeNode *> brickList = getAncestors() ;
  mat currentV = computeCovMat(m_knotsCoor, brickList.at(level)->GetKnotsCoor(), covParas) ;
  if (level > 0) {
    for (uint i = 0; i < level ; i++) {
      currentV -= m_bKnots.at(i) * brickList.at(i)->GetKmatrix() * trans(brickList.at(level)->GetBknots().at(i)) ;
    }
  }
  m_bKnots.at(level) = currentV ;
  if (GetChildren().at(0) != NULL) {
    for (auto & i : GetChildren()) {
      i->completeBknots(covParas, level) ;
    }
  }
}
