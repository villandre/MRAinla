// #ifdef _OPENMP
#include <omp.h>
// #endif

#include <math.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"
// #include "gperftools/profiler.h"

using namespace Rcpp ;
using namespace MRAinla ;
using namespace Eigen ;
using namespace std ;

struct gridPair{
  AugTree * grid ;
  vec vector ;
  gridPair() { } ;
  gridPair(AugTree * gridArg, vec vectorArg) : grid(gridArg), vector(vectorArg) { } ;
};

struct MVN{
  vec mu ;
  sp_mat precision ;

  MVN() { } ;
  MVN(const vec & mu, const sp_mat & precision): mu(mu), precision(precision) {} ;

  double ComputeExpTerm(const vec & coordinate) {
    vec output = -0.5 * (coordinate - mu).transpose() * precision * (coordinate - mu) ;
    return output(0) ;
  }
};

AugTree::AugTree(uint & Mlon,
                 uint & Mlat,
                 const spaceDimensions & mapDimensions,
                 inputdata & dataset,
                 inputdata & predictData,
                 const unsigned int numKnotsRes0,
                 double J,
                 const string & distMethod,
                 std::mt19937_64 & randomNumGenerator)
  : m_distMethod(distMethod), m_mapDimensions(mapDimensions), m_Mlon(Mlon), m_Mlat(Mlat)
{
  m_M = Mlon + Mlat ;

  m_assignedPredToKnot = Array<bool, Dynamic, 1>(predictData.spatialCoords.rows()) ;
  m_assignedPredToKnot.segment(0, m_assignedPredToKnot.size()) = false ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

  for (auto & i : tipNodes) {
    i->SetPredictLocations(predictData) ;
  }

  BuildTree(numKnotsRes0, J, dataset, predictData, randomNumGenerator) ;
  Rcout << "Grid built!" << std::endl ;

  m_numKnots = 0 ;

  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    m_numKnots += m_vertexVector.at(i)->GetKnotsCoor().spatialCoords.rows() ;
  }

  m_numObs = 0 ;
  for (auto & node : m_vertexVector) {
    m_numObs += node->GetObsInNode().size() ;
  }

}

void AugTree::BuildTree(const unsigned int numKnots0, double J, const inputdata & dataset, const inputdata & predictData, std::mt19937_64 & randomNumGenerator)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, dataset) ;

  m_vertexVector.push_back(topNode) ;
  ArrayXi numSplits(2) ;
  numSplits << m_Mlon, m_Mlat ;
  createLevels(topNode, "longitude", numSplits, dataset) ;
  numberNodes() ;

  generateKnots(topNode, numKnots0, J, randomNumGenerator, dataset, predictData) ;
}

void AugTree::numberNodes() {
  m_vertexVector.at(0)->SetNodeId(0) ;
  int currentIndex = 1 ;
  for (uint i = 1; i <= m_M; i++) {
    std::vector<TreeNode *> levelNodes = GetLevelNodes(i) ;
    for (auto & j : levelNodes) {
      j->SetNodeId(currentIndex) ;
      currentIndex += 1 ;
    }
  }
}

// We make sure that splits don't result in empty regions

void AugTree::createLevels(TreeNode * parent, std::string splitWhat, ArrayXi numSplitsLeft, const inputdata & dataset) {
  ArrayXi obsForMedian = parent->GetObsInNode() ;
  ArrayXi childMembership = ArrayXi::Zero(obsForMedian.size()) ;
  std::vector<spaceDimensions> childDimensions ;
  childDimensions.push_back(parent->GetDimensions()) ;

  if (obsForMedian.size() <= 1) {
    throw Rcpp::exception("Cannot split empty region or region with only one observation.\n") ;
  }

  // Handling longitude
  if (splitWhat == "longitude") {
    ArrayXi elementsInChild = ArrayXi::LinSpaced(obsForMedian.size(), 0, obsForMedian.size()-1) ;
    ArrayXd column = dataset.spatialCoords.col(0) ;
    ArrayXd elementsForMedian = elem(column, obsForMedian) ;
    double colMedian = median(elementsForMedian) ;
    ArrayXd updatedLongitude(2) ;
    updatedLongitude(0) = childDimensions.at(0).longitude(0) ;
    updatedLongitude(1) = colMedian ;
    ArrayXd newChildLongitude(2) ;
    newChildLongitude(0) = colMedian ;
    newChildLongitude(1) = childDimensions.at(0).longitude(1) ;
    spaceDimensions newDimensions = childDimensions.at(0) ;
    newDimensions.longitude = newChildLongitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(0).longitude = updatedLongitude ;
  } else if (splitWhat == "latitude") {
    ArrayXi elementsInChild = ArrayXi::LinSpaced(obsForMedian.size(), 0, obsForMedian.size()-1) ;
    ArrayXd column = dataset.spatialCoords.col(1) ;
    ArrayXd elementsForMedian = elem(column, obsForMedian) ;
    double colMedian = median(elementsForMedian) ;
    ArrayXd updatedLatitude(2) ;
    updatedLatitude(0) = childDimensions.at(0).latitude(0) ;
    updatedLatitude(1) = colMedian ;
    ArrayXd newChildLatitude(2) ;
    newChildLatitude(0) = colMedian ;
    newChildLatitude(1) = childDimensions.at(0).latitude(1) ;
    spaceDimensions newDimensions = childDimensions.at(0) ;
    newDimensions.latitude = newChildLatitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(0).latitude = updatedLatitude ;
  }

  int incrementedDepth = parent->GetDepth() + 1 ;
  for (auto & i : childDimensions) {
    if (parent->GetDepth() < m_M - 1) {
      InternalNode * newNode = new InternalNode(i, incrementedDepth, parent, dataset) ;
      m_vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    } else {
      TipNode * newNode = new TipNode(i, incrementedDepth, parent, dataset) ;
      m_numTips = m_numTips+1 ;
      m_vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    }

  }
  // Domain is split with respect to longitude, then latitude
  if (incrementedDepth < m_M) {
    if (splitWhat == "longitude") {
      if (numSplitsLeft(1) > 0) {
        splitWhat = "latitude" ;
        numSplitsLeft(1) -= 1 ;
      } else {
        numSplitsLeft(0) -= 1 ;
      }
    } else if (splitWhat == "latitude") {
      if (numSplitsLeft(0) > 0) {
        splitWhat = "longitude" ;
        numSplitsLeft(0) -= 1 ;
      } else {
        numSplitsLeft(1) -= 1 ;
      }
    }
    for (auto && i : parent->GetChildren()) {
      createLevels(i, splitWhat, numSplitsLeft, dataset) ;
    }
  }
}

void AugTree::generateKnots(TreeNode * node,
                            const unsigned int numKnotsRes0,
                            double J,
                            std::mt19937_64 & randomNumGenerator,
                            const inputdata & dataset,
                            const inputdata & predictData) {

  int numNodesAtLevel = GetLevelNodes(node->GetDepth()).size() ;
  int numKnotsToGen = std::max(uint(std::ceil((numKnotsRes0 * pow(J, node->GetDepth()))/numNodesAtLevel)), uint(2)) ;
  // node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetDepth() < m_M) {
    node->genKnotsOnSquare(predictData, numKnotsToGen, randomNumGenerator, m_assignedPredToKnot) ;
  } else {
    node->genKnotsOnSquare(dataset, numKnotsToGen, randomNumGenerator, m_assignedPredToKnot) ;
  }
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i, numKnotsRes0, J, randomNumGenerator, dataset, predictData) ;
    }
  }
}

void AugTree::computeWmats(const maternVec & MRAcovParasSpace, const double & spacetimeScaling, const double & spaceNuggetSD, const std::string & distMethod) {
  m_vertexVector.at(0)->ComputeWmat(MRAcovParasSpace, spacetimeScaling, spaceNuggetSD, distMethod) ;
  for (uint level = 1; level <= m_M; level++) {
    std::vector<TreeNode *> levelNodes = GetLevelNodes(level) ;

    // Trying openmp. We need to have a standard looping structure.
    // pragma omp parallel for

    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++)
    {
      (*it)->ComputeWmat(MRAcovParasSpace, spacetimeScaling, spaceNuggetSD, distMethod) ;
    }
  }
}

std::vector<TreeNode *> AugTree::GetLevelNodes(const uint & level) {
  std::vector<TreeNode *> nodesInLevel ;
  for (auto & i : m_vertexVector) {
    if (i->GetDepth() == level) {
      nodesInLevel.push_back(i) ;
    }
  }
  return nodesInLevel ;
}

std::vector<TreeNode *> AugTree::GetLevel(const uint level) {
  std::vector<TreeNode *> nodesAtLevel;
  if (level == 0) {
    nodesAtLevel.push_back(m_vertexVector.at(0)) ;
    return nodesAtLevel ;
  }
  for (auto & i : m_vertexVector) {
    if (i->GetDepth() == level) {
      nodesAtLevel.push_back(i) ;
    }
  }
  return nodesAtLevel ;
}

std::vector<TreeNode *> AugTree::Descendants(std::vector<TreeNode *> nodeList) {
  std::vector<TreeNode *> descendantList ;
  for (auto & i : nodeList) diveAndUpdate(i, &descendantList) ;
  return descendantList;
}

void AugTree::diveAndUpdate(TreeNode * nodePointer, std::vector<TreeNode *> * descendantList) {
  std::vector<TreeNode *> childrenVector = nodePointer->GetChildren() ;
  if (childrenVector.at(0) == NULL) {
    descendantList->push_back(nodePointer) ;
  } else {
    for (auto & i : childrenVector) diveAndUpdate(i, descendantList) ;
  }
}
