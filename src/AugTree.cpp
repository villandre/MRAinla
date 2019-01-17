#include <math.h>
#include <omp.h>
#include <gsl/gsl_randist.h>

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
  m_MRAcovParas.resize(1) ;
  m_MRAcovParas.fill(-1.5e6) ;

  gsl_rng_set(m_randomNumGenerator, seed) ;
  m_fixedEffParameters.resize(m_dataset.covariateValues.n_cols + 1) ; // An extra 1 is added to take into account the intercept.

  BuildTree(minObsForTimeSplit) ;
  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    m_vertexVector.at(i)->SetNodeId(i) ;
  }
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

double AugTree::ComputeMRAlogLik(const vec & thetaValues, const vec & MRAcovParas)
{
  cout << "Entering ComputeLoglik \n" ;
  bool sameCovParametersTest = false ;
  if (MRAcovParas.size() == m_MRAcovParas.size()) {
    sameCovParametersTest = approx_equal(m_MRAcovParas, MRAcovParas, "abs_diff", 1e-06) ;
  }
  printf("Equality test %i \n", sameCovParametersTest) ;
  if (!sameCovParametersTest) {
    m_MRAcovParas = MRAcovParas ;
    cout << "Entering computeWmats \n" ;
    m_MRAcovParas.print("MRA cov. parameters: ") ;
    computeWmats() ;
    cout << "Entering deriveAtildeMatrices \n" ;
    deriveAtildeMatrices() ;
  }
  cout << "Entering computeOmegas \n" ;
  computeOmegas(thetaValues) ;
  cout << "Entering computeU \n" ;
  computeU(thetaValues) ;
  cout << "Entering computeD \n" ;
  if (!sameCovParametersTest) {
    computeD() ;
  }
  cout << "Finalising evaluation \n" ;
  // The log-likelihood is a function of d and u computed at the root node, being at the
  // head of m_vertexVector, since it's the first one we ever created.
  double tempLogLik = (m_vertexVector.at(0)->GetD() + m_vertexVector.at(0)->GetU()) ;
  tempLogLik = -tempLogLik/2 ;
  return tempLogLik ;
}

void AugTree::computeWmats() {
  m_vertexVector.at(0)->ComputeBaseKmat(m_MRAcovParas) ;
  m_vertexVector.at(0)->ComputeWmat(m_MRAcovParas) ;

  for (uint level = 1; level < (m_M+1); level++) {
    std::vector<TreeNode *> levelNodes = getLevelNodes(level) ;
    // for (auto & i : levelNodes) {
    //   i->ComputeWmat() ;
    // }
    // Trying openmp. We need to have a standard looping structure.
    #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++)
    {
      (*it)->ComputeWmat(m_MRAcovParas) ;
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

void AugTree::computeOmegas(const vec & responseValues) {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveOmega(responseValues) ;
    }
  }
}

void AugTree::computeU(const vec & responseValues) {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveU(responseValues) ;
    }
  }
};

void AugTree::computeD() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;

    // Trying openmp. We need to have a standard looping structure.

    #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
      (*it)->DeriveD() ;
    }
  }
}

// The code only allows for main effects for now.

mat AugTree::ComputePosteriors(spatialcoor & predictionLocations, double & stepSize) {
  mat posteriorMatrix(predictionLocations.timeCoords.size() + m_dataset.covariateValues.n_cols + 1, 4, fill::zeros) ;

}

std::vector<GaussDistParas> AugTree::ComputeConditionalPrediction(const inputdata & predictionData) {
  distributePredictionData(predictionData) ;

  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = GetLevel(levelRecast) ;
    #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
    // for (auto & i : levelNodes) {
      (*it)->ComputeParasEtaDeltaTilde(predictionData, m_dataset, m_MRAcovParas) ;
    }
  }
  computeBtildeInTips() ;

  std::vector<vec> predictionsFromEachSection ;
  std::vector<TreeNode *> tipNodes = GetLevel(m_M) ;

  std::vector<GaussDistParas> distParasFromEachZone ;
  #pragma omp parallel for
  for (std::vector<TreeNode *>::iterator it = tipNodes.begin(); it < tipNodes.end(); it++) {
  // for (auto & i : tipNodes) {
    distParasFromEachZone.push_back((*it)->CombineEtaDelta(predictionData, m_fixedEffParameters)) ;
  }

  return distParasFromEachZone ;
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
      i->initiateBknots(m_MRAcovParas) ;
    }
  }
  for (uint level = 1; level < m_M ; level++) {
    for (auto & i : nodesAtLevels.at(level+1)) {
      i->completeBknots(m_MRAcovParas, level) ;
    }
  }
  for (auto & i : nodesAtLevels.at(m_M)) {
    i->computeBpred(m_predictLocations, m_MRAcovParas) ;
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
  vec intercept = ones<vec>(m_dataset.covariateValues.n_rows) ;
  mat incrementedCovar = m_dataset.covariateValues ;
  incrementedCovar.insert_cols(0, intercept) ;
  vec meanVec = vec(m_dataset.covariateValues.n_rows) ;
  meanVec = incrementedCovar * m_fixedEffParameters ;
  m_dataset.responseValues -= meanVec ;
}

arma::sp_mat AugTree::createSigmaStarInverse() {
  std::vector<mat *> KinverseAndBetaMatrixList = getKmatricesInversePointers() ;
  mat SigmaBetaInverse = eye<mat>(m_fixedEffParameters.size(), m_fixedEffParameters.size())/m_MRAcovParas.at(0) ; // Make sure those parameters are named.
  KinverseAndBetaMatrixList.push_back(&SigmaBetaInverse) ;
  sp_mat SigmaStarInverse = createSparseMatrix(KinverseAndBetaMatrixList) ;
  return SigmaStarInverse ;
}

arma::sp_mat AugTree::createHstar() {
  sp_mat Fmatrix = createFmatrix() ;
  vec intercept = ones<vec>(m_dataset.covariateValues.n_rows) ;
  mat covarCopy = m_dataset.covariateValues ;
  covarCopy.insert_cols(0, intercept) ;
  sp_mat incrementedCovar = conv_to<sp_mat>::from(covarCopy) ;

  sp_mat Hstar = join_rows(Fmatrix, incrementedCovar) ;
  return Hstar ;
}

arma::sp_mat AugTree::createFmatrix() {
  uint numKnotsTotal = 0 ;
  cout << "Creating F matrix... \n" ;
  for (auto & i : m_vertexVector) {
    numKnotsTotal += i->GetKnotsCoor().timeCoords.size() ;
  }
  sp_mat Fmat(m_dataset.timeCoords.size(), numKnotsTotal) ;
  std::vector<TreeNode *> tipNodes = getLevelNodes(m_M) ;

  uvec verticalBoundaryVec(m_M, fill::zeros) ;

  for (uint resolution = 1 ; resolution <= m_M ; resolution++) {
    std::vector<TreeNode *> levelNodes = getLevelNodes(m_M) ;
    uint numKnotsAtLevel = 0 ;
    for (auto & nodes : levelNodes) {
      numKnotsAtLevel += nodes->GetKnotsCoor().timeCoords.size() ;
    }
    verticalBoundaryVec(resolution) = verticalBoundaryVec(resolution - 1) + numKnotsAtLevel;
  }

  for (uint resolution = 0; resolution <= m_M ; resolution++) {
    for (auto & i : tipNodes) {

    }
  }
  cout << "Returning F matrix... \n" ;
  return Fmat ;
}

// For now, we assume that all hyperpriors have an inverse gamma distribution with the same parameters.

double AugTree::ComputeLogPriors(const arma::vec & MRAhyperPriorVec, const double errorSD, const double fixedParamSD, const double hyperAlpha, const double hyperBeta) {
  vec hyperPriorVec = MRAhyperPriorVec ;
  hyperPriorVec.resize(hyperPriorVec.size() + 2) ;
  hyperPriorVec(MRAhyperPriorVec.size()) = errorSD ;
  hyperPriorVec(MRAhyperPriorVec.size() + 1) = fixedParamSD ;
  vec invertedValues = 1/hyperPriorVec ;
  vec logContributions(invertedValues.size(), fill::zeros) ;

  std::transform(invertedValues.begin(), invertedValues.end(), logContributions.begin(),
                 [hyperAlpha, hyperBeta] (double & invertedValue) {return log(gsl_ran_gamma_pdf(invertedValue, hyperAlpha, hyperBeta)) + 2*log(invertedValue) ;}) ;
  return sum(logContributions) ;
}

double AugTree::ComputeLogJointCondTheta(const arma::vec & MRAvalues, const arma::vec & MRAcovParameters, const arma::vec & fixedEffParams, const double fixedEffSD) {

  double MRAlogLik = ComputeMRAlogLik(MRAvalues, MRAcovParameters) ;
  double fixedEffMean = 0 ;
  double fixedEffLogLik = 0 ;
  for (auto & i : fixedEffParams) {
    fixedEffLogLik += log(normpdf(i, fixedEffMean, fixedEffSD)) ;
  }
  return (MRAlogLik + fixedEffLogLik) ;
}


double AugTree::ComputeLogFullConditional(const arma::vec & MRAvalues, const arma::vec & fixedEffCoefs) {
  cout << "Creating Hstar... \n" ;
  sp_mat Hstar = createHstar() ;
  cout << "Creating SigmaStarInverse... \n" ;
  sp_mat SigmaStarInverse = createSigmaStarInverse() ;
  cout << "Creating Tmatrix... \n" ;
  sp_mat TmatrixInverse(m_dataset.timeCoords.size(), m_dataset.timeCoords.size()) ;
  TmatrixInverse.diag().fill(1/m_MRAcovParas.at(0)) ;
  cout << "Creating Qmat... \n" ;
  sp_mat Qmat = trans(Hstar) * TmatrixInverse * Hstar + SigmaStarInverse ;
  cout << "Computing bVec... \n" ;
  vec bVec = trans(Hstar) * TmatrixInverse * m_dataset.responseValues + SigmaStarInverse ;
  mat QmatInverse = inv_sympd(conv_to<mat>::from(Qmat)) ;
  vec meanVec = QmatInverse * bVec ;
  vec thetaValues = join_rows(MRAvalues, fixedEffCoefs) ;
  vec centeredThetas = thetaValues - meanVec ;
  double val;
  double sign;
  log_det(val, sign, QmatInverse) ;
  cout << "Finalising density evaluation..." ;
  double logDensPart1 = - 0.5 * QmatInverse.n_rows * (log((double) 2) +
                          log(M_PI)) + 0.5 * sign * val ;
  mat logDensPart2 = 0.5 * trans(centeredThetas) * QmatInverse * centeredThetas ;
  return logDensPart1 - logDensPart2(0,0) ;
}

double AugTree::ComputeGlobalLogLik(const arma::vec & MRAvalues, const arma::vec & fixedEffParams, const double errorSD) {
  mat incrementedCovar = m_dataset.covariateValues ;
  incrementedCovar.insert_cols(0, 1) ;
  incrementedCovar.col(0).fill(1) ;
  vec meanVec(m_dataset.responseValues.size(), fill::zeros) ;
  meanVec = incrementedCovar * fixedEffParams + MRAvalues ;
  vec sdVec(m_dataset.responseValues.size(), fill::zeros) ;
  sdVec.fill(errorSD) ;
  vec logDensVec(m_dataset.responseValues.size(), fill::zeros) ;
  logDensVec = log(normpdf(m_dataset.responseValues, meanVec, sdVec)) ;
  return sum(logDensVec) ;
}

double AugTree::ComputeJointPsiMarginal(const arma::vec & MRAhyperparas, const double errorSD, const double fixedEffSD, const double hyperAlpha, const double hyperBeta) {
  // In theory, this function should not depend on the theta values...
  // We can therefore arbitrarily set them all to 0.
  vec MRAvalues(m_dataset.responseValues.size(), fill::zeros) ;
  vec fixedEffCoefs(m_dataset.responseValues.size(), fill::zeros) ;
  double totalLogLik = ComputeGlobalLogLik(MRAvalues, fixedEffCoefs, errorSD) ;
  double logPrior = ComputeLogPriors(MRAhyperparas, errorSD, fixedEffSD, hyperAlpha, hyperBeta) ;
  double logCondDist = ComputeLogJointCondTheta(MRAvalues, MRAhyperparas, fixedEffCoefs, fixedEffSD) ;
  double logFullCond = ComputeLogFullConditional(MRAvalues, fixedEffCoefs) ;
  return logPrior + logCondDist + totalLogLik - logFullCond ;
}
