#include <gsl/gsl_sf_exp.h>
#include <math.h>

#include "TreeNode.h"

using namespace Eigen ;
using namespace MRAinla ;
using namespace std ;

ArrayXi TreeNode::deriveObsInNode(const spatialcoor & dataset) {
  Eigen::Array<bool, Dynamic, 1> lonCheck = (dataset.spatialCoords.col(0) > m_dimensions.longitude.matrix().minCoeff()) *
    (dataset.spatialCoords.col(0) <= m_dimensions.longitude.matrix().maxCoeff()) ; // Longitude check
  Eigen::Array<bool, Dynamic, 1> latCheck = (dataset.spatialCoords.col(1) > m_dimensions.latitude.matrix().minCoeff()) *
    (dataset.spatialCoords.col(1) <= m_dimensions.latitude.matrix().maxCoeff()) ; // Latitude check
  Eigen::Array<bool, Dynamic, 1> timeCheck = (dataset.timeCoords > m_dimensions.time.matrix().minCoeff()) *
    (dataset.timeCoords <= m_dimensions.time.matrix().maxCoeff()) ; // Time check
  ArrayXi output = find(lonCheck * latCheck * timeCheck) ;
  return output ;
}

// We assume a squared exponential decay function with sigma^2 = 1

double TreeNode::SqExpCovFunction(const Spatiotemprange & distance, const double & rhoSpace, const double & rhoTime, const double & spaceNuggetSD, const double & timeNuggetSD) {
  double spExp, timeExp ;

  if (distance.sp == 0) {
    spExp = 1 + spaceNuggetSD ;
  } else {
    spExp = exp(-pow(distance.sp, 2)/(2 * pow(rhoSpace, 2))) ;
  }

  if (distance.time == 0) {
    timeExp = 1 + timeNuggetSD ;
  } else {
    timeExp = exp(-pow(distance.time, 2)/(2 * pow(rhoTime, 2))) ;
  }
  return spExp * timeExp ;
}

double TreeNode::MaternCovFunction(const Spatiotemprange & distance, const maternVec & covParasSpace,
                                   const maternVec & covParasTime, const double & scaling, const double & spaceNuggetSD,
                                   const double & timeNuggetSD) {
  // Matern covariance with smoothness 1 and scaling 1.
  double spExp = maternCov(distance.sp, covParasSpace.m_rho, covParasSpace.m_smoothness , 1,
                           spaceNuggetSD) ;
  double timeExp = maternCov(distance.time, covParasTime.m_rho, covParasTime.m_smoothness, 1,
                             timeNuggetSD) ;
  return pow(scaling, 2) * spExp * timeExp ;
}

void TreeNode::baseComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {

  std::vector<TreeNode *> brickList = getAncestors() ;

  m_Wlist.at(0) = computeCovMat(m_knotsCoor, brickList.at(0)->GetKnotsCoor(), covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;

  for (uint l = 1; l <= m_depth; l++) {
    mat firstMat = computeCovMat(m_knotsCoor, brickList.at(l)->GetKnotsCoor(), covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;
    mat secondMat = mat::Zero(firstMat.rows(), firstMat.cols()) ;
    for (uint k = 0; k < l ; k++) {
      if (k < m_depth) {
      secondMat.noalias() += m_Wlist.at(k) *
        brickList.at(k)->GetKmatrix() *
        brickList.at(l)->GetWlist().at(k).transpose() ;
      } else {
        secondMat.noalias() += m_Wlist.at(m_depth) *
          brickList.at(k)->GetKmatrix() *
          brickList.at(l)->GetWlist().at(k).transpose() ;
      }
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

mat TreeNode::computeCovMat(const spatialcoor & spTime1, const spatialcoor & spTime2, const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {

  mat covMat = mat::Zero(spTime1.timeCoords.size(), spTime2.timeCoords.size()) ;
  for (uint rowIndex = 0; rowIndex < spTime1.timeCoords.size() ; rowIndex++) {
    for (uint colIndex = 0; colIndex < spTime2.timeCoords.size() ; colIndex++) {
      ArrayXd space1 = spTime1.spatialCoords.row(rowIndex) ;
      double time1 = spTime1.timeCoords(rowIndex) ;
      ArrayXd space2 = spTime2.spatialCoords.row(colIndex) ;
      double time2 = spTime2.timeCoords(colIndex) ;

      Spatiotemprange rangeValue = sptimeDistance(space1, time1, space2, time2) ;
      if (matern) {
        covMat(rowIndex, colIndex) = MaternCovFunction(rangeValue, covParasSp, covParasTime, scaling, spaceNuggetSD, timeNuggetSD) ;
      } else {
        covMat(rowIndex, colIndex) = SqExpCovFunction(rangeValue, covParasSp.m_rho, covParasTime.m_rho, spaceNuggetSD, timeNuggetSD) ;
      }
    }
  }

  return covMat ;
}
