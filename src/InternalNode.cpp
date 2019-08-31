#include "InternalNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void InternalNode::RemoveChild(TreeNode * childToRemove)
{
  auto ChildIterPos = std::find(m_children.begin(), m_children.end(), childToRemove) ;
  if (ChildIterPos == m_children.end())
  {
    cerr << "Warning: Trying to remove a child that was not found! \n" ;
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

    uvec aConverted(numKnots, fill::zeros) ;

    gsl_ran_choose(RNG, a, numKnots, b, dataCoor.timeCoords.size(), sizeof (double));

    for (uint i = 0 ; i < numKnots ; i++) {
      aConverted.at(i) = a[i] ;
    }

    mat jitteredSpace = dataCoor.spatialCoords.rows(aConverted) ;
    vec jitteredTime = dataCoor.timeCoords.elem(aConverted) ;

    for (auto & i : jitteredSpace) {
      i  += gsl_ran_gaussian(RNG, 0.001) ;
    }

    for(auto & i : jitteredTime) {
      i += gsl_ran_gaussian(RNG, 0.001) ;
    }
  } else {
    mat knotsSp(numKnots, 2, fill::zeros) ;

    float minLon = min(m_dimensions.longitude) ;
    float maxLon = max(m_dimensions.longitude) ;

    float minLat = min(m_dimensions.latitude) ;
    float maxLat = max(m_dimensions.latitude) ;

    float minTime = min(m_dimensions.time) ;
    float maxTime = max(m_dimensions.time) ;

    vec time(numKnots) ;

    knotsSp(0, 0) = (minLon + maxLon)/2 ;
    knotsSp(0, 1) = (minLat + maxLat)/2 ;
    time(0) = (minTime + maxTime)/2 ;

    // uint numCubes = ceil(double(numKnots - 1)/8) ;
    uint cubeRadiusInPoints = ceil(double(pow(numKnots, 1/3))) ;

    // double shortestLonRadius = (maxLon - minLon)/(numCubes + 1) ;
    // double shortestLatRadius = (maxLat - minLat)/(numCubes + 1) ;
    // double shortestTimeRadius = (maxTime - minTime)/(numCubes + 1) ;
    double lonDist = (maxLon - minLon)/cubeRadiusInPoints ;
    double latDist = (maxLat - minLat)/cubeRadiusInPoints ;
    double timeDist = (maxTime - minTime)/cubeRadiusInPoints ;

    uint rowIndex = 1 ; // The first knot in each zone is in the exact center of the zone.

    for (uint lonIndex = 0 ; lonIndex < cubeRadiusInPoints ; lonIndex++) {
      for (uint latIndex = 0 ; latIndex < cubeRadiusInPoints ; latIndex++) {
        for (uint timeIndex = 0 ; latIndex < cubeRadiusInPoints ; latIndex++) {
          knotsSp(rowIndex, 0) = minLon + double(lonIndex) * lonDist + gsl_ran_gaussian(RNG, 0.001) ;
          knotsSp(rowIndex, 1) = minLat + double(latIndex) * latDist + gsl_ran_gaussian(RNG, 0.001) ;
          time(rowIndex) = minTime + double(timeIndex) * timeDist  + gsl_ran_gaussian(RNG, 0.001) ;
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
      mat containerMat(m_children.at(0)->GetAtildeList(k, l).n_rows,
                           m_children.at(0)->GetAtildeList(k, l).n_cols, fill::zeros) ;

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

  m_Ktilde = inv_sympd(m_KtildeInverse) ;

  for (uint k = 0; k <= m_depth ; k++) {
    for (uint l = 0; l <= k ; l++) {

      m_AtildeList.at(k).at(l) = m_Alist.at(k).at(l) -
        trans(m_Alist.at(m_depth).at(k)) * m_Ktilde * m_Alist.at(m_depth).at(l) ;
    }
  }
}

void InternalNode::DeriveOmega(const arma::vec & responseValues) {
  vec containerVec ;
  for (uint k = 0; k <= m_depth ; k++) {
    containerVec.set_size(m_children.at(0)->GetOmegaTilde(k).size()) ;
    containerVec.fill(0) ;

    for (auto & i: m_children) {
      containerVec += i->GetOmegaTilde(k) ;
    }
    m_omega.at(k) = containerVec ;
  }

  // for (auto & i : m_children) {
  //   i->clearOmegaTilde() ;
  // }

  for (uint k = 0; k <= m_depth ; k++) {
    vec secondTerm = trans(m_Alist.at(m_depth).at(k)) * m_Ktilde * m_omega.at(m_depth) ;
    m_omegaTilde.at(k) = m_omega.at(k) -  secondTerm;
  }
}

void InternalNode::DeriveU(const arma::vec & responseValues) {
  mat firstTerm = -trans(m_omega.at(m_depth)) * m_Ktilde * m_omega.at(m_depth) ;
  double secondTerm = 0 ;
  // secondTerm = std::accumulate(m_children.begin(), m_children.end(), secondTerm,
  //                              [](double a, TreeNode * b) {
  //                                return a + b->GetU();
  //                              }) ;
  for (auto & i : m_children) {
    secondTerm += i->GetU() ;
  }
  m_u = firstTerm.at(0,0) + secondTerm ;
}

void::InternalNode::DeriveD() {
  double value = 0 ;
  double sign = 0 ;
  log_det(value, sign, m_KtildeInverse) ; // Function is supposed to update value and sign in place.
  // sign is supposed to be positive, since a negative determinant would indicate negative variance!
  if (sign < 0) {
    throw Rcpp::exception("Precision matrix KtildeInverse cannot have negative determinant. \n") ;
  }
  m_d = value ;
  log_det(value, sign, GetKmatrixInverse()) ;
  if (sign < 0) {
    throw Rcpp::exception("Precision matrix Kinverse cannot have negative determinant. \n") ;
  }
  m_d = m_d - value ;
  double thirdTerm = 0 ;
  // thirdTerm = std::accumulate(m_children.begin(), m_children.end(), thirdTerm,
  //                              [](double a, TreeNode * b) {
  //                                return a + b->GetD();
  //                              }) ;
  for (auto & i : m_children) {
    thirdTerm += i->GetD() ;
  }
  m_d += thirdTerm ;
}

void InternalNode::ComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {
  baseComputeWmat(covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;
  m_Wlist.at(m_depth) = symmatl(m_Wlist.at(m_depth)) ;
  m_K = inv_sympd(GetKmatrixInverse()) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
}

// void InternalNode::ComputeParasEtaDeltaTilde(const spatialcoor & predictLocations, const inputdata & dataset, const arma::vec & covParas) {
//   m_etaTilde.meanPara = m_Ktilde * m_omega.at(m_depth) ;
//   m_etaTilde.covPara = m_Ktilde ; // This is repetition. It should not harm memory too much though.
// }

// void InternalNode::ComputeBaseKmat(const vec & covParaVec) {
//   mat covMatrix = ComputeCovMatrix(covParaVec) ;
//   m_Kinverse = symmatl(covMatrix) ;
//   m_K = inv_sympd(covMatrix) ;
// }
