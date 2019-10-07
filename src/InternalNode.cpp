#include "InternalNode.h"
#include <gsl/gsl_randist.h>

using namespace Eigen ;
using namespace MRAinla ;

void InternalNode::RemoveChild(TreeNode * childToRemove)
{
  auto ChildIterPos = std::find(m_children.begin(), m_children.end(), childToRemove) ;
  if (ChildIterPos == m_children.end())
  {
    Rcpp::warning("Warning: Trying to remove a child that was not found! \n") ;
  }
  else
  {
    m_children.erase(ChildIterPos) ;
  }
}

void InternalNode::genRandomKnots(spatialcoor & dataCoor, const uint & numKnots, const gsl_rng * RNG) {

  if (m_depth == 0) {
    double a[numKnots], b[dataCoor.timeCoords.size()] ;

    for (uint i = 0; i < dataCoor.timeCoords.size() ; i++)
    {
      b[i] = (double) i;
    }

    ArrayXi aConverted = ArrayXi::Zero(numKnots) ;

    gsl_ran_choose(RNG, a, numKnots, b, dataCoor.timeCoords.size(), sizeof (double));

    for (uint i = 0 ; i < numKnots ; i++) {
      aConverted(i) = a[i] ;
    }

    ArrayXXd jitteredSpace = rows(dataCoor.spatialCoords, aConverted) ;
    ArrayXd jitteredTime = elem(dataCoor.timeCoords, aConverted) ;

    for (uint i = 0; i < jitteredSpace.size(); i++) {
      jitteredSpace(i) += gsl_ran_gaussian(RNG, 0.0001) ;
    }

    for (uint i = 0; i < jitteredTime.size(); i++) {
      jitteredTime(i) += gsl_ran_gaussian(RNG, 0.0001) ;
    }
    m_knotsCoor = spatialcoor(jitteredSpace, jitteredTime) ;
  } else {
    ArrayXXd knotsSp = ArrayXXd::Zero(numKnots, 2) ;

    double minLon = m_dimensions.longitude.minCoeff() ;
    double maxLon = m_dimensions.longitude.maxCoeff() ;

    double minLat = m_dimensions.latitude.minCoeff() ;
    double maxLat = m_dimensions.latitude.maxCoeff() ;

    double minTime = m_dimensions.time.minCoeff() ;
    double maxTime = m_dimensions.time.maxCoeff() ;

    ArrayXd time = ArrayXd::Zero(numKnots) ;

    double cubeRadiusInPoints = ceil(double(cbrt(numKnots))) ;

    double offsetPerc = 0.01 ;
    double lonDist = (maxLon - minLon) * (1-offsetPerc * 2)/(cubeRadiusInPoints - 1) ;
    double latDist = (maxLat - minLat) * (1-offsetPerc * 2)/(cubeRadiusInPoints - 1) ;
    double timeDist = (maxTime - minTime) * (1-offsetPerc * 2)/(cubeRadiusInPoints - 1) ;

    uint rowIndex = 0 ;
    for (uint lonIndex = 0 ; lonIndex < cubeRadiusInPoints ; lonIndex++) {
      for (uint latIndex = 0 ; latIndex < cubeRadiusInPoints ; latIndex++) {
        for (uint timeIndex = 0 ; timeIndex < cubeRadiusInPoints ; timeIndex++) {
          knotsSp(rowIndex, 0) = minLon + offsetPerc * (maxLon - minLon) + double(lonIndex) * lonDist + gsl_ran_gaussian(RNG, 0.0001) ;
          knotsSp(rowIndex, 1) = minLat + offsetPerc * (maxLat - minLat) + double(latIndex) * latDist + gsl_ran_gaussian(RNG, 0.0001) ;
          time(rowIndex) = minTime + offsetPerc * (maxTime - minTime) + double(timeIndex) * timeDist  + gsl_ran_gaussian(RNG, 0.0001) ;
          rowIndex += 1 ;
          if (rowIndex >= numKnots) break ;
        }
        if (rowIndex >= numKnots) break ;
      }
      if (rowIndex >= numKnots) break ;
    }
    m_knotsCoor = spatialcoor(knotsSp, time) ;
  }
}

void InternalNode::DeriveAtilde() {

  for (uint k = 0; k <= m_depth ; k++) {
    for (uint l = 0; l <= k ; l++) {
      mat containerMat = mat::Zero(m_children.at(0)->GetAtildeList(k, l).rows(),
                           m_children.at(0)->GetAtildeList(k, l).cols()) ;

      for (auto & i : m_children) {
        containerMat += i->GetAtildeList(k,l) ;
      }
      m_Alist.at(k).at(l) = containerMat ;
    }
  }
  // for (auto & i : m_children) {
  //   i->clearAtildeList() ;
  // }
  m_KtildeInverse = GetKmatrixInverse() + m_Alist.at(m_depth).at(m_depth) ;

  m_Ktilde = m_KtildeInverse.ldlt().solve(mat::Identity(m_KtildeInverse.rows(), m_KtildeInverse.rows())) ;

  for (uint k = 0; k <= m_depth ; k++) {
    for (uint l = 0; l <= k ; l++) {

      m_AtildeList.at(k).at(l) = m_Alist.at(k).at(l) -
        m_Alist.at(m_depth).at(k).transpose() * m_Ktilde * m_Alist.at(m_depth).at(l) ;
    }
  }
}

void InternalNode::DeriveOmega(const vec & responseValues) {
  for (uint k = 0; k <= m_depth ; k++) {
    vec containerVec = vec::Zero(m_children.at(0)->GetOmegaTilde(k).size()) ;

    for (auto & i: m_children) {
      containerVec += i->GetOmegaTilde(k) ;
    }
    m_omega.at(k) = containerVec ;
  }

  for (uint k = 0; k <= m_depth ; k++) {
    vec secondTerm = m_Alist.at(m_depth).at(k).transpose() * m_Ktilde * m_omega.at(m_depth) ;
    m_omegaTilde.at(k) = m_omega.at(k) -  secondTerm;
  }
}

void InternalNode::DeriveU(const vec & responseValues) {
  mat firstTerm = -m_omega.at(m_depth).transpose() * m_Ktilde * m_omega.at(m_depth) ;
  double secondTerm = 0 ;

  for (auto & i : m_children) {
    secondTerm += i->GetU() ;
  }
  m_u = firstTerm(0,0) + secondTerm ;
}

void::InternalNode::DeriveD() {

  double value = m_KtildeInverse.ldlt().vectorD().array().log().sum() ;

  m_d = value ;
  double otherValue = GetKmatrixInverse().ldlt().vectorD().array().log().sum() ;
  m_d = m_d - otherValue ;
  double thirdTerm = 0 ;

  for (auto & i : m_children) {
    thirdTerm += i->GetD() ;
  }
  m_d += thirdTerm ;
}

void InternalNode::ComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {
  baseComputeWmat(covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;

  // m_Wlist.at(m_depth).triangularView<Upper>() = m_Wlist.at(m_depth).triangularView<Lower>() ; // Will this cause aliasing?

  m_K = GetKmatrixInverse().ldlt().solve(mat::Identity(GetKmatrixInverse().rows(), GetKmatrixInverse().cols())) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
}
