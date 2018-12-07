#include <cmath>
#include <omp.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"

using namespace Rcpp ;
using namespace arma ;
using namespace MRAinla ;

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, vec & observations, mat & obsSp, uvec & obsTime, uint & minObsForTimeSplit, unsigned long int & seed, mat & covariates)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime, covariates) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;

  gsl_rng_set(m_randomNumGenerator, seed) ;
  m_fixedEffParameters.resize(m_dataset.covariateValues.n_cols + 1) ; // An extra 1 is added to take into account the intercept.

  BuildTree(minObsForTimeSplit) ;
}

void AugTree::BuildTree(const uint & minObsForTimeSplit)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, m_dataset) ;

  m_vertexVector.push_back(topNode) ;

  createLevels(topNode, minObsForTimeSplit) ;

  generateKnots(topNode) ;
}

void AugTree::createLevels(TreeNode * parent, const uint & numObsForTimeSplit) {
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
      m_numTips = m_numTips+1 ;
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

void AugTree::generateKnots(TreeNode * node) {
  // float scaledRespSize =  0.05 ;
  // float depthContribution = 19*node->GetDepth() ;
  // depthContribution = depthContribution/60 ;
  // float expToRound = depthContribution + scaledRespSize ; // Why is it that we have to store every term separately for expToRound to not evaluate to 0?
  // float numKnots = expToRound * node->GetObsInNode().size() ;
  //// CAREFUL: THIS IS HARD-CODED, THERE MIGHT BE A BETTER WAY /////////////
  float numKnots = sqrt((float) node->GetObsInNode().size()) ;
  uint numKnotsToGen = uint (std::ceil(numKnots)) ;
  node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i) ;
    }
  }
}

void AugTree::ComputeMRAloglik()
{
  cout << "Entering ComputeLoglik \n" ;
  computeWmats() ;
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
  double tempLogLik = (m_vertexVector.at(0)->GetD() + m_vertexVector.at(0)->GetU()) ;
  tempLogLik = -tempLogLik/2 ;
  m_MRAlogLik = tempLogLik ;
}

void AugTree::computeWmats() {
  m_vertexVector.at(0)->ComputeBaseKmat(m_covParameters) ;
  m_vertexVector.at(0)->ComputeWmat(m_covParameters) ;

  for (uint level = 1; level < (m_M+1); level++) {
    std::vector<TreeNode *> levelNodes = getLevelNodes(level) ;
    // for (auto & i : levelNodes) {
    //   i->ComputeWmat() ;
    // }
    // Trying openmp. We need to have a standard looping structure.
    #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++)
    {
      (*it)->ComputeWmat(m_covParameters) ;
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

void AugTree::deriveAtildeMatrices() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;

    // for (auto & i : levelNodes) {
    //   i->DeriveAtilde() ;
    // }
    #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
      (*it)->DeriveAtilde() ;
    }
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

    // Trying openmp. We need to have a standard looping structure.

    // #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
      (*it)->DeriveD() ;
    }
  }
}

// The code only allows for main effects for now.

mat AugTree::ComputePosteriors(spatialcoor & predictionLocations, double & stepSize) {
  mat posteriorMatrix(predictionLocations.timeCoords.size() + m_dataset.covariateValues.n_cols + 1, 4, fill::zeros) ;

}

std::vector<GaussDistParas> AugTree::ComputeConditionalPrediction(const spatialcoor & predictionLocations) {
  distributePredictionData(predictionLocations) ;

  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = GetLevel(levelRecast) ;
    for (auto & i : levelNodes) {
      i->ComputeParasEtaDeltaTilde(predictionLocations, m_dataset, m_covParameters) ;
    }
  }
  computeBtildeInTips() ;

  std::vector<vec> predictionsFromEachSection ;
  std::vector<TreeNode *> tipNodes = GetLevel(m_M) ;
  std::vector<GaussDistParas> distParasFromEachZone(tipNodes.size()) ;
  // std::transform(tipNodes.begin(), tipNodes.end(), distParasFromEachZone, [] (TreeNode * node) {
  //   return node->CombineEtaDelta() ; // Check the ordering of observations.
  // }) ;

  std::vector<GaussDistParas>::iterator outputIterator = distParasFromEachZone.begin() ;
  for (auto & i : tipNodes) {
    *outputIterator = i->CombineEtaDelta(m_dataset, m_fixedEffParameters) ;
    std::next(outputIterator) ;
  }
  return distParasFromEachZone ;
}

double AugTree::ComputeGlobalLogLik() {
  // First we create the necessary GSL vectors...
  mat incrementedCovar(m_dataset.covariateValues) ;
  incrementedCovar.insert_cols(0, 1) ;
  incrementedCovar.col(0).fill(1) ;
  vec meanVec(m_dataset.responseValues.size(), 0) ;
  meanVec = incrementedCovar * m_fixedEffParameters + m_spatialComponents ;
  vec sdVec(m_dataset.responseValues.size(), m_errorSD) ;
  vec logDensVec(m_dataset.responseValues.size(), 0) ;
  logDensVec = log(normpdf(m_dataset.responseValues, meanVec, sdVec)) ;
  return sum(logDensVec) ;
}

void AugTree::distributePredictionData(const spatialcoor & predictLocations) {
  m_predictLocations = predictLocations ;

  for (auto & i : m_vertexVector) {
    i->SetPredictLocations(predictLocations) ;
  }
}

void AugTree::computeBtildeInTips() {
  std::vector<std::vector<TreeNode *>> nodesAtLevels(m_M+1) ;
  for (uint i = 0 ; i <= m_M ; i++) {
    nodesAtLevels.at(i) = GetLevel(i) ;
  }
  // The two following loops cannot be joined.
  for (uint level = 1; level <= m_M ; level++) {
    for (auto & i : nodesAtLevels.at(level)) {
      i->initiateBknots(m_covParameters) ;
    }
  }

  for (uint level = 1; level < m_M ; level++) {
    for (auto & i : nodesAtLevels.at(level+1)) {
      i->completeBknots(m_covParameters, level) ;
    }
  }

  for (auto & i : nodesAtLevels.at(m_M)) {
    i->computeBpred(m_predictLocations, m_covParameters) ;
    i->deriveBtilde(m_predictLocations) ;
  }
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

void AugTree::CleanPredictionComponents() {
  //TO_DO
}

void AugTree::CenterResponse() {
  cout << "Entering CenterResponse() \n" ;
  mat incrementedCovar = m_dataset.covariateValues ;
  incrementedCovar.insert_cols(0, 1, TRUE) ;
  incrementedCovar.col(0).fill(1) ;
  vec meanVec(m_dataset.responseValues.size(), 0) ;
  meanVec = incrementedCovar * m_fixedEffParameters ;
  m_dataset.responseValues -= meanVec ;
  cout << "Leaving CenterResponse() \n" ;
}
