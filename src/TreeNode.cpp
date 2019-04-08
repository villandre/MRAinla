#include <gsl/gsl_sf_exp.h>
#include <math.h>

#include "TreeNode.h"

using namespace arma ;
using namespace MRAinla ;
using namespace std ;

arma::uvec TreeNode::deriveObsInNode(const spatialcoor & dataset) {
  uvec lonCheck = (dataset.spatialCoords.col(0) > min(m_dimensions.longitude)) %
    (dataset.spatialCoords.col(0) <= max(m_dimensions.longitude)) ; // Longitude check
  uvec latCheck = (dataset.spatialCoords.col(1) > min(m_dimensions.latitude)) %
    (dataset.spatialCoords.col(1) <= max(m_dimensions.latitude)) ; // Latitude check
  uvec timeCheck = (dataset.timeCoords > min(m_dimensions.time)) %
    (dataset.timeCoords <= max(m_dimensions.time)) ; // Time check
  return find(lonCheck % latCheck % timeCheck) ; // find is equivalent to which in R
}

// We assume a squared exponential decay function with sigma^2 = 1

double TreeNode::SqExpCovFunction(const Spatiotemprange & distance, const vec & covParameters, const double & spaceNuggetSD, const double & timeNuggetSD) {
  double spExp, timeExp ;
  if (distance.sp == 0) {
    spExp = pow(spaceNuggetSD, 2) ;
  } else {
    spExp = exp(-pow(distance.sp, 2)/(2 * pow(covParameters.at(0), 2))) ;
  }

  if (distance.time == 0) {
    timeExp = pow(timeNuggetSD, 2) ;
  } else {
    timeExp = exp(-pow(distance.time, 2)/(2 * pow(covParameters.at(1), 2))) ;
  }
  return spExp * timeExp ;
}

double TreeNode::MaternCovFunction(const Spatiotemprange & distance, const vec & covParameters, const double & spaceNuggetSD, const double & timeNuggetSD) {
  // Matern covariance with smoothness 1 and scaling 1.
  double spExp = maternCov(distance.sp, covParameters.at(0), 1 , 1, spaceNuggetSD) ;
  double timeExp = maternCov(distance.sp, covParameters.at(1), 1 , 1, timeNuggetSD) ;
  return spExp * timeExp ;
}

void TreeNode::baseComputeWmat(const vec & covParas, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {
  std::vector<TreeNode *> brickList = getAncestors() ;

  m_Wlist.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor(), covParas, matern, spaceNuggetSD, timeNuggetSD) ;

  for (uint l = 1; l <= m_depth; l++) {
    mat firstMat = computeCovMat(m_knotsCoor, brickList.at(l)->GetKnotsCoor(), covParas, matern, spaceNuggetSD, timeNuggetSD) ;
    mat secondMat(firstMat.n_rows, firstMat.n_cols, fill::zeros) ;
    for (uint k = 0; k < l ; k++) {
      secondMat += m_Wlist.at(k) *
        brickList.at(k)->GetKmatrix() *
        trans(brickList.at(l)->GetWlist().at(k)) ;
    }
    m_Wlist.at(l) = firstMat - secondMat ;
    if (l == m_depth) {
      double sign1, sign2 ;
      double value1, value2 ;
      log_det(value1, sign1, firstMat) ;
      log_det(value2, sign2, secondMat) ;
    }
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

mat TreeNode::computeCovMat(const spatialcoor & spTime1, const spatialcoor & spTime2, const vec & covParas, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {

  mat covMat(spTime1.timeCoords.size(), spTime2.timeCoords.size(), fill::zeros) ;
  for (uint rowIndex = 0; rowIndex < spTime1.timeCoords.size() ; rowIndex++) {
    for (uint colIndex = 0; colIndex < spTime2.timeCoords.size() ; colIndex++) {
      vec space1 = conv_to<vec>::from(spTime1.spatialCoords.row(rowIndex)) ;
      double time1 = spTime1.timeCoords.at(rowIndex) ;
      vec space2 = conv_to<vec>::from(spTime2.spatialCoords.row(colIndex)) ;
      double time2 = spTime2.timeCoords.at(colIndex) ;

      Spatiotemprange rangeValue = sptimeDistance(space1, time1, space2, time2) ;
      if (matern) {
        covMat.at(rowIndex, colIndex) = MaternCovFunction(rangeValue, covParas, spaceNuggetSD, timeNuggetSD) ;
      } else {
        covMat.at(rowIndex, colIndex) = SqExpCovFunction(rangeValue, covParas, spaceNuggetSD, timeNuggetSD) ;
      }
    }
  }

  return covMat ;
}

// void TreeNode::initiateBknots(const vec & covParas) {
//   std::vector<TreeNode *> brickList = getAncestors() ;
//   if (m_depth > 0) {
//     m_bKnots.resize(m_depth) ;
//     m_bKnots.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor(), covParas) ;
//   }
// }

// void TreeNode::completeBknots(const vec & covParas, const uint level) {
//   std::vector<TreeNode *> brickList = getAncestors() ;
//   mat currentV = computeCovMat(m_knotsCoor, brickList.at(level)->GetKnotsCoor(), covParas) ;
//   if (level > 0) {
//     for (uint i = 0; i < level ; i++) {
//       currentV -= m_bKnots.at(i) * brickList.at(i)->GetKmatrix() * trans(brickList.at(level)->GetBknots().at(i)) ;
//     }
//   }
//   m_bKnots.at(level) = currentV ;
//   if (GetChildren().at(0) != NULL) {
//     for (auto & i : GetChildren()) {
//       i->completeBknots(covParas, level) ;
//     }
//   }
// }
