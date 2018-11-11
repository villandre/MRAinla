#include <math.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"
#include "TreeNode.h"

using namespace Rcpp ;
using namespace arma ;
using namespace MRAinla ;

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, vec & observations, mat & obsSp, uvec & obsTime, uint & minObsForTimeSplit, unsigned long int & seed)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  gsl_rng_set(m_randomNumGenerator, seed) ;

  BuildTree(minObsForTimeSplit) ;
}

void AugTree::BuildTree(uint & minObsForTimeSplit)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, m_dataset) ;
  m_vertexVector.push_back(topNode) ;

  createLevels(topNode, minObsForTimeSplit) ;
  cout << "We cleared createLevels! \n" ;
  generateKnots(topNode) ;
}

void AugTree::createLevels(TreeNode * parent, uint & numObsForTimeSplit) {
  uvec obsForMedian(m_dataset.responseValues.size(), fill::zeros) ;
  obsForMedian = parent->GetObsInNode() ;

  if (obsForMedian.n_elem == 0) {
    throw Rcpp::exception("Empty region.") ;
  }
  mat spMediansMat = median(m_dataset.spatialCoords.rows(obsForMedian), 0) ;
  rowvec spMedians = spMediansMat.row(0) ;

  uvec sortedTimes = sort(m_dataset.timeCoords.elem(obsForMedian)) ;
  uint timeMedian = sortedTimes.at(std::ceil((obsForMedian.size()-1)/2));

  std::vector<dimensions> dimensionsForChildren ;
  vec lonRange(2, fill::zeros), latRange(2, fill::zeros) ;

  if (parent->GetObsInNode().size() < numObsForTimeSplit) {
    for (uint i = 0; i < 2; i++) {
      for (uint j = 0 ; j < 2; j++) {
        lonRange.at(0) = parent->GetDimensions().longitude.at(i) ;
        lonRange.at(1) = spMedians.at(0);
        lonRange = sort(lonRange);

        latRange.at(0) = parent->GetDimensions().latitude.at(j) ;
        latRange.at(1) = spMedians.at(1);
        latRange = sort(latRange);

        dimensionsForChildren.push_back(dimensions(lonRange, latRange, parent->GetDimensions().time)) ;
      }
    }
  } else {
    uvec timeRange(2, fill::zeros) ;
    for (uint i = 0; i < 2; i++) {
      for (uint j = 0; j < 2; j++) {
        for(uint k = 0; k < 2; k++) {
          lonRange.at(0) = parent->GetDimensions().longitude.at(i) ;
          lonRange.at(1) = spMedians.at(0);
          lonRange = sort(lonRange);

          latRange.at(0) = parent->GetDimensions().latitude.at(j) ;
          latRange.at(1) = spMedians.at(1);
          latRange = sort(latRange);

          timeRange.at(0) = parent->GetDimensions().time.at(k) ;
          timeRange.at(1) = timeMedian ;
          timeRange = sort(timeRange) ;

          dimensionsForChildren.push_back(dimensions(lonRange, latRange, timeRange)) ;
        }
      }
    }
  }
  uint incrementedDepth = parent->GetDepth()+1 ;

  for (auto &&i : dimensionsForChildren) {
    if (parent->GetDepth() < m_M-1) {
      InternalNode * newNode = new InternalNode(i, incrementedDepth, parent, m_dataset) ;
      m_vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    } else {
      TipNode * newNode = new TipNode(i, incrementedDepth, parent, m_dataset) ;
      m_vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    }
  }
  if (incrementedDepth < m_M) {
    incrementedDepth = incrementedDepth + 1 ;
    for (auto && i : parent->GetChildren()) {
      createLevels(i, numObsForTimeSplit) ;
    }
  }
}

void AugTree::ComputeLoglik(const std::vector<mat> & withinClusTransProbs, const std::vector<mat> & betweenClusTransProbs, const vec & limProbs)
{
// TO_DO
}

void AugTree::generateKnots(TreeNode * node) {
  cout << "Entering generate knots... \n \n" ;
  double expToRound = (node->GetDepth()*19)/60 + 1/20 * m_dataset.responseValues.size() ; //// PROBLEME HERE!
  cout << "expToRound: " << expToRound << "\n " ;
  uint numKnotsToGen = std::ceil(expToRound) ;
  node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i) ;
    }
  }
  cout << "Leaving generateKnots... \n \n" ;
}
