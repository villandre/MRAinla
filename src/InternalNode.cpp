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

    m_knotsCoor = spatialcoor(spaceCoords) ;
  }

  numKnots = min(numKnots, int(m_obsInNode.size())) ;

  m_knotsCoor = spatialcoor(ArrayXXd::Zero(numKnots, 2)) ;

  double minLon = m_dimensions.longitude.minCoeff() ;
  double maxLon = m_dimensions.longitude.maxCoeff() ;

  double minLat = m_dimensions.latitude.minCoeff() ;
  double maxLat = m_dimensions.latitude.maxCoeff() ;

  std::uniform_real_distribution<double> distribution(0, 1);
  for (uint i = 0; i < m_knotsCoor.spatialCoords.rows(); i++) {
    m_knotsCoor.spatialCoords(i,0) =  distribution(generator) * (maxLon - minLon) + minLon ;
    m_knotsCoor.spatialCoords(i,1) =  distribution(generator) * (maxLat - minLat) + minLat ;
  }
}

void InternalNode::genKnotsOnSquare(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & generator, Array<bool, Dynamic, 1> & assignedPredLocations) {
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

  if (numUnassignedPreds > 0) {
    numKnotsFromPred = min(numKnots, numUnassignedPreds) ;

    std::shuffle(unassignedPredObsInNode.begin(), unassignedPredObsInNode.end(), generator) ;
    ArrayXi shuffledPredIndices(numKnotsFromPred) ;
    for (uint i = 0 ; i < numKnotsFromPred; i++) {
      shuffledPredIndices(i) = unassignedPredObsInNode.at(i) ;
      assignedPredLocations(unassignedPredObsInNode.at(i)) = true ;
    }

    spacePredCoords = rows(dataCoor.spatialCoords, shuffledPredIndices) ;
  }

  ArrayXXd knotsSp(0, 2) ;

  if (numKnots - numKnotsFromPred > 0) {
    numKnots -= numKnotsFromPred ;

    int numPointsPerEdge = ceil((double(numKnots) + 16)/12) ;
    // For now, let's check what happens when the number of knots is always 8 + 12*i,
    // thus ensuring a uniform spread in the space.
    numKnots = numPointsPerEdge * 12 - 16 ;
    knotsSp.resize(numKnots, 2) ;

    knotsSp = ArrayXXd::Zero(numKnots, 2) ;

    double minLon = m_dimensions.longitude.minCoeff() ;
    double maxLon = m_dimensions.longitude.maxCoeff() ;

    double minLat = m_dimensions.latitude.minCoeff() ;
    double maxLat = m_dimensions.latitude.maxCoeff() ;

    double offsetPerc = 1e-5 ;
    double lonDist = (maxLon - minLon) * (1-offsetPerc * 2)/(double(numPointsPerEdge) - 1) ;
    double latDist = (maxLat - minLat) * (1-offsetPerc * 2)/(double(numPointsPerEdge) - 1) ;

    Array2d scalingArray = {lonDist, latDist} ;
    Array2d offsetCoords = {minLon + offsetPerc * (maxLon - minLon),
                            minLat + offsetPerc * (maxLat - minLat)
                            } ;

    intSquare squareForKnots(numPointsPerEdge, scalingArray, offsetCoords) ;

    // We place the corners first.

    std::array<Eigen::Array2d, 4> corners = squareForKnots.getCorners(false) ;
    uint index = 0 ;
    for (auto & corner : corners) {
      knotsSp.row(index) = corner.head(2) ;
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
        Array2d knotArray = squareForKnots.getEdgePointCoor(edgeIndex, pointIndex) ;
        knotsSp.row(index) = knotArray.head(2) ;
        edgeIndex += 1 ;
        index += 1 ;
      }
    }
  }
  ArrayXXd mergedSpace = join_cols(spacePredCoords, knotsSp) ;

  std::uniform_real_distribution<double> distribution(-1e-6, 1e-6);
  for (uint i = 0; i < mergedSpace.rows(); i++) {
    mergedSpace(i) += distribution(generator) ;
    mergedSpace(i + mergedSpace.rows()) += distribution(generator) ;
  }
  m_knotsCoor = spatialcoor(mergedSpace) ;
}

void InternalNode::ComputeWmat(const maternVec & covParasSp, const double & scaling, const double & spaceNuggetSD, const string & distMethod) {
  baseComputeWmat(covParasSp, scaling, spaceNuggetSD, distMethod) ;

  // m_Wlist.at(m_depth).triangularView<Upper>() = m_Wlist.at(m_depth).triangularView<Lower>() ; // Will this cause aliasing?

  m_K = GetKmatrixInverse().selfadjointView<Upper>().ldlt().solve(mat::Identity(GetKmatrixInverse().rows(), GetKmatrixInverse().cols())) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
}
