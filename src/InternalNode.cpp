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

void InternalNode::genRandomKnots(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & generator) {

  if (m_depth == 0) {
    ArrayXXd spaceCoords(0, 2) ;
    ArrayXd timeCoords(0) ;

    std::vector<int> obsInNodeVec ;

    for (uint i = 0; i < m_obsInNode.size(); i++) {
      obsInNodeVec.push_back(m_obsInNode(i)) ;
    }

    std::shuffle(obsInNodeVec.begin(), obsInNodeVec.end(), generator) ;

    ArrayXi shuffledIndices(numKnots) ;
    for (uint i = 0 ; i < numKnots; i++) {
      shuffledIndices(i) = obsInNodeVec.at(i) ;
    }
    spaceCoords = rows(dataCoor.spatialCoords, shuffledIndices) ;
    timeCoords = elem(dataCoor.timeCoords, shuffledIndices) ;

    m_knotsCoor = spatialcoor(spaceCoords, timeCoords) ;
  }

  numKnots = min(numKnots, int(m_obsInNode.size())) ;

  m_knotsCoor = spatialcoor(ArrayXXd::Zero(numKnots, 2), ArrayXd::Zero(numKnots)) ;

  double minLon = m_dimensions.longitude.minCoeff() ;
  double maxLon = m_dimensions.longitude.maxCoeff() ;

  double minLat = m_dimensions.latitude.minCoeff() ;
  double maxLat = m_dimensions.latitude.maxCoeff() ;

  double minTime = m_dimensions.time.minCoeff() ;
  double maxTime = m_dimensions.time.maxCoeff() ;

  std::uniform_real_distribution<double> distribution(0, 1);
  for (uint i = 0; i < m_knotsCoor.spatialCoords.rows(); i++) {
    m_knotsCoor.spatialCoords(i,0) =  distribution(generator) * (maxLon - minLon) + minLon ;
    m_knotsCoor.spatialCoords(i,1) =  distribution(generator) * (maxLat - minLat) + minLat ;
    m_knotsCoor.timeCoords(i) =  distribution(generator) * (maxTime - minTime) + minTime ;
  }
}


void InternalNode::genKnotsOnCube(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & generator, Array<bool, Dynamic, 1> & assignedPredLocations) {
  numKnots = min(numKnots, int(m_obsInNode.size())) ;
  ArrayXi predObsInNode = deriveObsInNode(dataCoor) ;

  std::vector<int> unassignedPredObsInNode ;
  if (predObsInNode.size() > 0) {
    for (uint i = 0; i < predObsInNode.size(); i++) {
      if (!assignedPredLocations(predObsInNode(i)))
        unassignedPredObsInNode.push_back(predObsInNode(i)) ;
    }
  }

  int numUnassignedPreds = unassignedPredObsInNode.size() ;

  int numKnotsFromPred = 0 ;
  ArrayXXd spacePredCoords(0, 2) ;
  ArrayXd timePredCoords(0) ;

  if (numUnassignedPreds > 0) {
    numKnotsFromPred = min(numKnots, numUnassignedPreds) ;

    std::shuffle(unassignedPredObsInNode.begin(), unassignedPredObsInNode.end(), generator) ;
    ArrayXi shuffledPredIndices(numKnotsFromPred) ;
    for (uint i = 0 ; i < numKnotsFromPred; i++) {
      shuffledPredIndices(i) = unassignedPredObsInNode.at(i) ;
      assignedPredLocations(unassignedPredObsInNode.at(i)) = true ;
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

    double offsetPerc = 1e-5 ;
    double lonDist = (maxLon - minLon) * (1-offsetPerc * 2)/(double(numPointsPerEdge) - 1) ;
    double latDist = (maxLat - minLat) * (1-offsetPerc * 2)/(double(numPointsPerEdge) - 1) ;
    double timeDist = (maxTime - minTime) * (1-offsetPerc * 2)/(double(numPointsPerEdge) - 1) ;

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

  std::uniform_real_distribution<double> distribution(-1e-6, 1e-6);
  for (uint i = 0; i < mergedTime.size(); i++) {
    mergedTime(i) += distribution(generator) ;
    mergedSpace(i) += distribution(generator) ;
    mergedSpace(i + mergedSpace.rows()) += distribution(generator) ;
  }
  m_knotsCoor = spatialcoor(mergedSpace, mergedTime) ;
}

void InternalNode::ComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const double & nuggetSD, const string & distMethod) {
  baseComputeWmat(covParasSp, covParasTime, scaling, nuggetSD, distMethod) ;

  // m_Wlist.at(m_depth).triangularView<Upper>() = m_Wlist.at(m_depth).triangularView<Lower>() ; // Will this cause aliasing?

  m_K = GetKmatrixInverse().selfadjointView<Upper>().ldlt().solve(mat::Identity(GetKmatrixInverse().rows(), GetKmatrixInverse().cols())) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
}
