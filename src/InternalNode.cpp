#include "InternalNode.h"

using namespace Eigen ;
using namespace MRAinla ;
using namespace std ;

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

void InternalNode::genRandomKnots(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & generator, Array<bool, Dynamic, 1> & assignedPredLocations) {

  numKnots = min(numKnots, int(m_obsInNode.size())) ;
  ArrayXi predObsInNode = deriveObsInNode(dataCoor) ;
  std::vector<int> subAssignedPredLocations ;
  if (predObsInNode.size() > 0) {
    for (uint i = 0; i < predObsInNode.size(); i++) {
      subAssignedPredLocations.push_back(assignedPredLocations(predObsInNode(i))) ;
    }
  }

  int numUnassignedPreds = 0 ;

  for (auto & i : subAssignedPredLocations) {
    if (!i) numUnassignedPreds += 1 ;
  }
  int numKnotsFromPred = 0 ;
  ArrayXXd spacePredCoords(0, 2) ;
  ArrayXd timePredCoords(0) ;
  if (numUnassignedPreds > 0) {
    numKnotsFromPred = min(numKnots, numUnassignedPreds) ;

    std::vector<int> vectorToShuffle ;
    for (uint i = 0; i < subAssignedPredLocations.size(); i++) {
      if (!subAssignedPredLocations.at(i)) {
        vectorToShuffle.push_back(i) ;
      }
    }
    std::shuffle(vectorToShuffle.begin(), vectorToShuffle.end(), generator) ;
    ArrayXi shuffledPredIndices(numKnotsFromPred) ;
    for (uint i = 0 ; i < numKnotsFromPred; i++) {
      shuffledPredIndices(i) = vectorToShuffle.at(i) ;
      assignedPredLocations(shuffledPredIndices(i)) = true ;
    }

    spacePredCoords = rows(dataCoor.spatialCoords, shuffledPredIndices) ;
    timePredCoords = elem(dataCoor.timeCoords, shuffledPredIndices) ;
  }

  ArrayXXd knotsSp(0, 2) ;
  ArrayXd time(0) ;

  if (numKnots - numKnotsFromPred > 0) {
    numKnots -= numKnotsFromPred ;

    int numPointsPerEdge = ceil((double(numKnots) + 16)/12) ;
    // For now, let's check what happens when the number of knots is always 8 + 12*i,
    // thus ensuring a uniform spread in the space.
    numKnots = numPointsPerEdge * 12 - 16 ;
    knotsSp.resize(numKnots, 2) ;
    time.resize(numKnots) ;

    knotsSp = ArrayXXd::Zero(numKnots, 2) ;

    double minLon = m_dimensions.longitude.minCoeff() ;
    double maxLon = m_dimensions.longitude.maxCoeff() ;

    double minLat = m_dimensions.latitude.minCoeff() ;
    double maxLat = m_dimensions.latitude.maxCoeff() ;

    double minTime = m_dimensions.time.minCoeff() ;
    double maxTime = m_dimensions.time.maxCoeff() ;

    double offsetPerc = 0.00001 ;
    double lonDist = (maxLon - minLon) * (1-offsetPerc * 2)/(numPointsPerEdge - 1) ;
    double latDist = (maxLat - minLat) * (1-offsetPerc * 2)/(numPointsPerEdge - 1) ;
    double timeDist = (maxTime - minTime) * (1-offsetPerc * 2)/(numPointsPerEdge - 1) ;

    Array3d scalingArray = {lonDist, latDist, timeDist} ;
    Array3d offsetCoords = {minLon + offsetPerc * (maxLon - minLon),
                            minLat + offsetPerc * (maxLat - minLat),
                            minTime + offsetPerc * (maxTime - minTime)} ;

    intCube cubeForKnots(numPointsPerEdge, scalingArray, offsetCoords) ;

    // We place the corners first.

    std::array<Eigen::Array3d, 8> corners = cubeForKnots.getCorners(false) ;
    uint index = 0 ;
    for (auto & corner : corners) {
      knotsSp.row(index) = corner.head(2) ;
      time(index) = corner(2) ;
      index += 1 ;
    }

    if (numKnots > 8) {
      uint numRemainingKnots = numKnots - 8 ;
      uint pointIndex = 0 ;
      uint edgeIndex = 0 ;
      for (uint i = 0 ; i < numRemainingKnots; i++) {
        if ((edgeIndex % 12) == 0) {
          edgeIndex = 0 ;
          pointIndex += 1 ;
        }
        Array3d knotArray = cubeForKnots.getEdgePointCoor(edgeIndex, pointIndex) ;
        knotsSp.row(index) = knotArray.head(2) ;
        time(index) = knotArray(2) ;
        edgeIndex += 1 ;
        index += 1 ;
      }
    }
  }
  ArrayXXd mergedSpace = join_cols(spacePredCoords, knotsSp) ;
  ArrayXd mergedTime = join_cols(timePredCoords, time) ;

  std::uniform_real_distribution<double> distribution(0.0, 0.000001);
  for (uint i = 0; i < mergedTime.size(); i++) {
    mergedTime(i) += distribution(generator) ;
    mergedSpace(i) += distribution(generator) ;
    mergedSpace(i + mergedTime.size()) += distribution(generator) ;
  }
  m_knotsCoor = spatialcoor(mergedSpace, mergedTime) ;
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
  double otherValue = GetKmatrixInverse().selfadjointView<Upper>().ldlt().vectorD().array().log().sum() ;
  m_d = m_d - otherValue ;
  double thirdTerm = 0 ;

  for (auto & i : m_children) {
    thirdTerm += i->GetD() ;
  }
  m_d += thirdTerm ;
}

void InternalNode::ComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const double & spaceNuggetSD, const double & timeNuggetSD, const string & distMethod) {
  baseComputeWmat(covParasSp, covParasTime, scaling, spaceNuggetSD, timeNuggetSD, distMethod) ;

  // m_Wlist.at(m_depth).triangularView<Upper>() = m_Wlist.at(m_depth).triangularView<Lower>() ; // Will this cause aliasing?

  m_K = GetKmatrixInverse().selfadjointView<Upper>().ldlt().solve(mat::Identity(GetKmatrixInverse().rows(), GetKmatrixInverse().cols())) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
}
