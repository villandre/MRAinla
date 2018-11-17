#include <math.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"

using namespace Rcpp ;
using namespace arma ;
using namespace MRAinla ;

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, vec & observations, mat & obsSp, uvec & obsTime, uint & minObsForTimeSplit, unsigned long int & seed, vec & covPars)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  gsl_rng_set(m_randomNumGenerator, seed) ;

  BuildTree(minObsForTimeSplit, covPars) ;
}

void AugTree::BuildTree(const uint & minObsForTimeSplit, const vec & covPara)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, m_dataset, covPara) ;
  m_vertexVector.push_back(topNode) ;

  createLevels(topNode, minObsForTimeSplit, covPara) ;

  generateKnots(topNode) ;
}

void AugTree::createLevels(TreeNode * parent, const uint & numObsForTimeSplit, const vec & covPara) {
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
  float scaledRespSize =  0.05 ;
  float depthContribution = 19*node->GetDepth() ;
  depthContribution = depthContribution/60 ;
  float expToRound = depthContribution + scaledRespSize ; // Why is it that we have to store every term separately for expToRound to not evaluate to 0?
  float numKnots = expToRound * node->GetObsInNode().size() ;

  uint numKnotsToGen = uint (std::ceil(numKnots)) ;
  node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i) ;
    }
  }
}

void AugTree::ComputeLoglik()
{
  cout << "Entering ComputeLoglik \n" ;
  computeWmats() ;
  cout << "Entering setBtips \n" ;
  setBtips() ;
  cout << "Entering setSigmaTips \n" ;
  setSigmaTips() ;
  cout << "Entering deriveAtildeMatrices \n" ;
  deriveAtildeMatrices() ;
  cout << "Entering computeOmegas \n" ;
  computeOmegas() ;
  cout << "Entering computeU \n" ;
  computeU() ;
  cout << "Entering computeD \n" ;
  computeD() ;
  cout << "Finalising evaluation \n" ;
  // The log-likelihood is a function of d and u computed at the root node, being at the
  // head of m_vertexVector, since it's the first one we ever created.
  m_logLik = (m_vertexVector.at(0)->GetD() + m_vertexVector.at(0)->GetU())/2 ;
}

void AugTree::computeWmats() {
  m_vertexVector.at(0)->ComputeBaseKmat() ;
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

void AugTree::setBtips() {
  std::vector<TreeNode *> allTips = getLevelNodes(m_M) ;
  for (auto && i : allTips) {
    i->DeriveB() ;
  }
}

void AugTree::deriveAtildeMatrices() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    printf("Deriving Atilde for level %u. \n", level) ;
    for (auto & i : levelNodes) {
      i->DeriveAtilde() ;
    }
    cout << "Done! \n" ;
  }
}

void AugTree::computeOmegas() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveOmega(m_dataset) ;
    }
  }
}

void AugTree::computeU() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveU(m_dataset) ;
    }
  }
};

void AugTree::computeD() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveD() ;
    }
  }
}
