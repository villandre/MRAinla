#include <math.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"
#include "TreeNode.h"

using namespace Rcpp ;
using namespace arma ;
using namespace MRAinla ;

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, vec & observations, mat & obsSp, uvec & obsTime, uint & minObsForTimeSplit, unsigned long int & seed, vec & covarianceParameter)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  gsl_rng_set(m_randomNumGenerator, seed) ;

  BuildTree(minObsForTimeSplit, covarianceParameter) ;
}

void AugTree::BuildTree(uint & minObsForTimeSplit, vec & covPara)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, m_dataset, covPara) ;
  m_vertexVector.push_back(topNode) ;

  createLevels(topNode, minObsForTimeSplit, covPara) ;
  cout << "We cleared createLevels! \n" ;
  generateKnots(topNode) ;
}

void AugTree::createLevels(TreeNode * parent, uint & numObsForTimeSplit, vec & covPara) {
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
      InternalNode * newNode = new InternalNode(i, incrementedDepth, parent, m_dataset, covPara) ;
      m_vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    } else {
      TipNode * newNode = new TipNode(i, incrementedDepth, parent, m_dataset, covPara) ;
      m_numTips = m_numTips+1 ;
      m_vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    }
  }
  if (incrementedDepth < m_M) {
    incrementedDepth = incrementedDepth + 1 ;
    for (auto && i : parent->GetChildren()) {
      createLevels(i, numObsForTimeSplit, covPara) ;
    }
  }
}

void AugTree::generateKnots(TreeNode * node) {
  cout << "Entering generate knots... \n \n" ;

  float scaledRespSize =  0.05 ;
  float depthContribution = 19*node->GetDepth() ;
  depthContribution = depthContribution/60 ;
  cout << "Depth contribution: " << depthContribution << "\n \n" ;
  float expToRound = depthContribution + scaledRespSize ; // Why is it that we have to store every term separately for expToRound to not evaluate to 0?
  float numKnots = expToRound * node->GetObsInNode().size() ;

  uint numKnotsToGen = uint (std::ceil(numKnots)) ;
  cout << "Number of knots: " << numKnotsToGen << "\n \n" ;
  node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i) ;
    }
  }
  cout << "Leaving generateKnots... \n \n" ;
}

void AugTree::ComputeLoglik()
{
  computeWmats() ;
  setBtips() ;
  setSigmaTips() ;
  setAtildeTips() ;
  recurseA() ;
  setOmegaTildeTips() ;
  recurseOmega() ;

  computeUtips() ;
  recurseU() ;

  computeDtips() ;
  recurseD() ;

  m_logLik <- (m_d + m_u)/2 ;
}

void AugTree::computeWmats() {
  mat baseK = m_vertexVector.at(0)->ComputeBaseKmat() ;
  m_vertexVector.at(0)->SetKandInv(baseK) ;
  m_vertexVector.at(0)->ComputeWmat() ;
  for (uint level = 1; level < (m_M+1); level++) {
    std::vector<TreeNode *> levelNodes = getLevelNodes(level) ;
    for (auto & i : levelNodes) {
      i->ComputeWmat() ;
    }
  }
}

std::vector<TreeNode *> AugTree::getLevelNodes(uint & level) {
  std::vector<TreeNode *> nodesInLevel ;
  for (auto & i : m_vertexVector) {
    if (i->GetDepth() == level) {
      nodesInLevel.push_back(i) ;
    }
  }
  return nodesInLevel ;
}

void AugTree::setSigmaTips() {

  std::vector<TreeNode *> allTips = getLevelNodes(m_M) ;

  for (auto & i : allTips) {
    i->SetSigma(i->GetWlist().at(m_M)) ;
  }
}

void AugTree::setAtildeTips() {
  std::vector<TreeNode *> allTips = getLevelNodes(m_M) ;

  for (auto & i : allTips) {
    i->DeriveAtilde() ;
  }
}

void AugTree::setBtips() {
  std::vector<TreeNode *> allTips = getLevelNodes(m_M) ;
  for (auto & i : allTips) {
    i->GetBlist().resize(m_M) ;
    auto startIter = i->GetWlist().begin() ;
    auto copyIter = i->GetBlist().begin() ;
    std::copy(startIter, startIter + m_M - 1, copyIter) ;
  }
}

void AugTree::recurseA() {
  for (int level = m_M - 1 ; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    mat containerMat ;
    for (auto & i : levelNodes) {
      for (uint k = 0; k <= i->GetDepth() ; k++) {
        for (uint l = 0; l <= k ; l++) {
          containerMat.reshape(i->GetChildren().at(0)->GetAtildeList().at(k).at(l).n_rows,
                               i->GetChildren().at(0)->GetAtildeList().at(k).at(l).n_cols) ;
          containerMat.fill(0) ;
          containerMat = std::accumulate(i->GetChildren().begin(), i->GetChildren().end(), containerMat,
                                         [](std::vector<std::vector<mat>> a, std::vector<std::vector<mat>> b) {
                                          return a.at(k).at(l) + b.at(k).at(l);
                                          }) ;

        }
      }

    }
  }
}
