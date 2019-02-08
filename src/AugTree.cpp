#include <math.h>
#include <omp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"

using namespace Rcpp ;
using namespace arma ;
using namespace MRAinla ;

struct gridPair{
  AugTree * grid ;
  vec * vector ;
};

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
  m_numKnots = 0 ;
  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    m_vertexVector.at(i)->SetNodeId(i) ;
    m_numKnots += m_vertexVector.at(i)->GetKnotsCoor().timeCoords.size() ;
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
    computeWmats() ;
    deriveAtildeMatrices() ;
  }
  computeOmegas(thetaValues) ;
  computeU(thetaValues) ;
  if (!sameCovParametersTest) {
    computeD() ;
  }
  // The log-likelihood is a function of d and u computed at the root node, being at the
  // head of m_vertexVector, since it's the first one we ever created.
  double tempLogLik = (m_vertexVector.at(0)->GetD() + m_vertexVector.at(0)->GetU()) ;
  tempLogLik = -tempLogLik/2 ;
  return tempLogLik ;
}

void AugTree::computeWmats() {
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
  sp_mat smallFmat = Fmatrix.cols(0, 100) ;

  vec intercept = ones<vec>(m_dataset.covariateValues.n_rows) ;
  mat covarCopy = m_dataset.covariateValues ;
  covarCopy.insert_cols(0, intercept) ;
  sp_mat incrementedCovar = conv_to<sp_mat>::from(covarCopy) ;
  sp_mat Hstar = join_rows(Fmatrix, incrementedCovar) ;
  return Hstar ;
}

// Very long fonction. Should probably be shortened...
// The F mat is a huge sparse matrix with number of rows (columns) equal to the total
// number of observations (knots).
// The function starts by inputting in Fmat matrices at resolution M (in a block-diagonal fashion).
// All matrices added to F are taken from m_BmatList(resolution) in the tip nodes,
// with resolution = M at this point.
// To make sure there's correspondence in rows between matrices inputted at all resolutions
// the function creates a parent matrix for all sets of siblings at the level currently considered.
// This matrix is a combination of m_BmatList(resolution - 1) across siblings.
// The matrices created are bound in a list.
// Once all sets of siblings at a given resolution have been processed, the algorithm moves
// up to resolution = resolution - 1.  The list of matrices created previously must now be processed.
// To do that, the algorithm follows the same steps outlined before (initially, the list
// considered was made up of m_BmatList(m_M) for all tip nodes).
// The algorithm ends when resolution 0 is reached.

arma::sp_mat AugTree::createFmatrix() {
  uint numKnotsTotal = 0 ;

  for (auto & i : m_vertexVector) {
    numKnotsTotal += i->GetNumKnots() ;
  }
  sp_mat Fmat(m_dataset.timeCoords.size(), numKnotsTotal) ;
  std::vector<TreeNode *> tipNodes = getLevelNodes(m_M) ;
  std::vector<std::pair<TreeNode *, mat>> currentLevelMatrices ;
  std::vector<std::pair<TreeNode *, mat>> parentLevelMatrices ;
  uint horizontalBoundary ;
  uint verticalBoundary = numKnotsTotal ; // Vertical boundary starts on the right of the very last knot.
  std::vector<TreeNode *> nodeSiblings ;
  std::pair<TreeNode *, mat> siblingParentPair ;

  for (auto & i : tipNodes) {
    currentLevelMatrices.push_back(std::make_pair(i, i->GetB(m_M))) ; // We start iterating at resolution M, hence the initialisation implemented here.
  }
  uint numProcessedKnots = 0 ;

  for (uint inverseResol = 0 ; inverseResol <= m_M ; inverseResol++) {
    uint resolution = m_M - inverseResol ;
    horizontalBoundary = 0 ;
    uint numKnotsAtLevel = 0 ;
    for (auto & j : GetLevel(resolution)) numKnotsAtLevel += j->GetNumKnots() ;
    numProcessedKnots += numKnotsAtLevel ;
    verticalBoundary = numKnotsTotal - numProcessedKnots ;

    while (currentLevelMatrices.size() > 0) { // currentLevelMatrices contains matrices that have not been added to Fmat yet. Once it is empty, it means they have all been added.
      if (resolution > 0) {
        nodeSiblings = currentLevelMatrices.at(0).first->Siblings() ;
        siblingParentPair = std::make_pair(nodeSiblings.at(0)->GetParent(),
          mat(0, nodeSiblings.at(0)->GetParent()->GetNumKnots(), fill::zeros)) ; // This will eventually contain information about the parent matrix for all the siblings.
      } else {
        nodeSiblings.push_back(currentLevelMatrices.at(0).first) ;
      }

      for (auto & i : nodeSiblings) {
        uint index = 0 ;
        bool check = currentLevelMatrices.at(0).first == i ;
        // Find the element of currentLevelMatrices that corresponds to the sibling currently being processed.
        while (!check) {
          check = currentLevelMatrices.at(index).first == i ;
          index += 1 ;
        }
        // Update Fmat
        Fmat(horizontalBoundary, verticalBoundary, size(currentLevelMatrices.at(index).second)) = currentLevelMatrices.at(index).second ;
        // Update placement indices for the blocks in Fmat.
        horizontalBoundary += i->GetObsInNode().size() ;
        verticalBoundary += i->GetNumKnots() ;
        // The node sibling is removed from currentLevelMatrices once it has been processed.
        currentLevelMatrices.erase(currentLevelMatrices.begin() + index) ;
      }
      // Prepare the list of matrices one resolution up
      if (resolution > 0) {
        std::vector<TreeNode *> descendedTips = Descendants(nodeSiblings) ;
        siblingParentPair.second = mat(siblingParentPair.first->GetNumKnots(), 0) ; // We'll transpose the result of binding matrices columnwise. This is faster than binding rowwise.
        for (auto & i : descendedTips) {
          // Matrices in Armadillo are column-major. Horizontal joins are therefore faster than vertical joins.
          siblingParentPair.second = join_horiz(siblingParentPair.second, trans(i->GetB(resolution - 1))) ; // Vertical join, like cbind in R, join_horiz would be faster, since matrix are column-based in Armadillo.
        }
        siblingParentPair.second = trans(siblingParentPair.second) ; // Transposing is fast.
        parentLevelMatrices.push_back(siblingParentPair) ;
        nodeSiblings.clear() ; // This vector will be re-incremented, hence the need to empty it.
      }
    }
    currentLevelMatrices = parentLevelMatrices ;
    parentLevelMatrices.clear() ; // parentLevelMatrices will be re-populated from scratch, hence the need to empty it.
  }
  return Fmat ;
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


double AugTree::ComputeLogFullConditional(const arma::vec & meanParaValues, const arma::vec & fixedEffCoefs) {

  sp_mat Hstar = createHstar() ;
  sp_mat SigmaStarInverse = createSigmaStarInverse() ;
  sp_mat TmatrixInverse(m_dataset.timeCoords.size(), m_dataset.timeCoords.size()) ;
  TmatrixInverse.diag().fill(1/m_MRAcovParas.at(0)) ;
  sp_mat Qmat = trans(Hstar) * TmatrixInverse * Hstar + SigmaStarInverse ;
  // The formulation for bVec is valid if priors for the eta's and fixed effect coefficients have mean zero, else, a second term comes into play Sigma * mu ;
  sp_mat bVec = trans(Hstar) * TmatrixInverse * conv_to<sp_mat>::from(m_dataset.responseValues) ;

  Eigen::SparseMatrix<double> eigen_s = Rcpp::as<Eigen::SparseMatrix<double>>(Rcpp::wrap(Qmat));
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(eigen_s);
  double ldet = solver.logAbsDeterminant();
  sp_mat thetaValues= sp_mat(meanParaValues) ;

  sp_mat logKernelResult = - 0.5 * (trans(thetaValues) * Qmat * thetaValues + trans(thetaValues) * bVec) ;

  double logFullValue = -0.5 * SigmaStarInverse.n_rows * log(2*PI) + 0.5 * ldet + logKernelResult(0,0) ;

  return logFullValue;
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

double AugTree::ComputeJointPsiMarginal(const arma::vec & MRAhyperparas, const double fixedEffSD, const double errorSD, const double hyperAlpha, const double hyperBeta) {
  // In theory, this function should not depend on the theta values...
  // We can therefore arbitrarily set them all to 0.
  cout << "Entering ComputeJointPsiMarginal... \n" ;
  vec MRAvalues(m_numKnots + m_dataset.covariateValues.n_cols + 1, fill::zeros) ;
  vec MRAvaluesInterpolant(m_dataset.timeCoords.n_rows, fill::zeros) ;
  vec fixedParValues(m_dataset.covariateValues.n_cols + 1, fill::zeros) ;
  cout << "Computing global log-lik... \n" ;
  double totalLogLik = ComputeGlobalLogLik(MRAvaluesInterpolant, fixedParValues, errorSD) ;
  cout << "Computing log-priors... \n" ;
  double logPrior = ComputeLogPriors(MRAhyperparas, errorSD, fixedEffSD, hyperAlpha, hyperBeta) ;
  cout << "Computing log-conditional... \n" ;
  double logCondDist = ComputeLogJointCondTheta(MRAvalues, MRAhyperparas, fixedParValues, fixedEffSD) ;
  cout << "Computing log-full conditional... \n" ;
  double logFullCond = ComputeLogFullConditional(MRAvalues, fixedParValues) ;
  cout << "Summing... \n" ;
  return logPrior + logCondDist + totalLogLik - logFullCond ;
}

// This inversion is based on recursive partitioning of the Q matrix. It is based on the observation that it is
// possible to form block-diagonal matrices on the diagonal which can be easily inverted.
// The challenge in inverting Qmat was the very heavy memory burden.
// This function involves much smaller matrices, which will make the operations easier to handle.
// mat AugTree::invertQmat(const sp_mat & Qmat) {
//   // We start in the lower-right corner and use M levels of nesting.
//   // We find all block diagonal matrices on the diagonal.
//   std::vector<uvec> blockDiagonalMatricesIndices ;
//   int lowerBound = Qmat.n_rows - 1 ;
//   cout << "Detecting block diagonal matrices... \n" ;
//   do {
//     uvec indicesInBlock = extractBlockIndicesFromLowerRight(Qmat.submat(0, 0, lowerBound, lowerBound)) ;
//     blockDiagonalMatricesIndices.push_back(indicesInBlock) ;
//     lowerBound = indicesInBlock.at(0) - 1 ;
//   } while (lowerBound >= 0) ;
//   cout << "Done! \n" ;
//   printf("Number of block diagonal matrices : %i \n", int(blockDiagonalMatricesIndices.size())) ;
//   for (auto & i : blockDiagonalMatricesIndices) {
//     printf("Number of blocks: %i \n" , int(i.size())) ;
//   }
//   cout << "Initializing D... \n" ;
//   uint numRowsD = blockDiagonalMatricesIndices.at(0).tail(1)(0) - blockDiagonalMatricesIndices.at(0)(0) ;
//   // printf("Indices first/last element in D: %i %i \n", blockDiagonalMatricesIndices.at(0)(0), blockDiagonalMatricesIndices.at(0).tail(1)(0)) ;
//   // blockDiagonalMatricesIndices.at(0).print("Block indices:") ;
//   uint shift = blockDiagonalMatricesIndices.at(0).at(0) ;
//   uvec shiftedBlockIndices(blockDiagonalMatricesIndices.size()) ;
//   shiftedBlockIndices = blockDiagonalMatricesIndices.at(0) - shift ;
//   cout << "Initializing Dinv... \n" ;
//   sp_mat Dinv = invertSymmBlockDiag(Qmat(blockDiagonalMatricesIndices.at(0).at(0),
//                                          blockDiagonalMatricesIndices.at(0).at(0),
//                                          size(numRowsD, numRowsD)),
//                                          shiftedBlockIndices) ;
//
//   // int reverseIndexA ;
//   cout << "Launching the recursion... \n" ;
//   for (uint index = 2 ; index < blockDiagonalMatricesIndices.size() ; index++) {
//     // reverseIndexA = blockDiagonalMatricesIndices.size() - index ;
//     uint nrowA =  blockDiagonalMatricesIndices.at(index).tail(1)(0) -
//       blockDiagonalMatricesIndices.at(index).at(0) ;
//     uint nrowD = Dinv.n_rows ;
//     printf("Number of rows of A and D: %i %i \n", int(nrowA), int(nrowD)) ;
//     printf("Processing A matrix %i \n", index) ;
//     shift = blockDiagonalMatricesIndices.at(index).at(0) ;
//     shiftedBlockIndices = blockDiagonalMatricesIndices.at(index) - shift ;
//     invFromDecomposition(Qmat(blockDiagonalMatricesIndices.at(index).at(0),
//                                        blockDiagonalMatricesIndices.at(index).at(0),
//                                        size(nrowA, nrowA)),
//                                   Qmat(blockDiagonalMatricesIndices.at(index).at(0),
//                                         blockDiagonalMatricesIndices.at(index).at(0) + nrowA,
//                                         size(nrowA, nrowD)),
//                                   Qmat(blockDiagonalMatricesIndices.at(index).at(0) + nrowA,
//                                         blockDiagonalMatricesIndices.at(index).at(0) + nrowA,
//                                         size(nrowD, nrowD)),
//                                         &Dinv, shiftedBlockIndices) ;
//     // D = Qmat(blockDiagonalMatricesIndices.at(index).at(0),
//     //            blockDiagonalMatricesIndices.at(index).at(0),
//     //            size(nrowA + nrowD, nrowA + nrowD)) ;
//   }
//   cout << "Recursion done! Exiting...\n" ;
//   return mat(1,1, fill::zeros) ; // PLACEHOLDER!!!
// }

// uvec AugTree::extractBlockIndicesFromLowerRight(const arma::sp_mat & symmSparseMatrix) {
//   std::vector<uint> blockIndices ;
//   blockIndices.push_back(symmSparseMatrix.n_rows) ; // We start by adding a bound on the right, although other points in the vector correspond to upper-left corners.
//   int index = symmSparseMatrix.n_rows - 2 ;
//   uint lowerBoundary = symmSparseMatrix.n_rows - 1 ;
//   uvec nonZeros ;
//   if (any(vec(symmSparseMatrix.diag()) == 0)) {
//     throw Rcpp::exception("Full conditional has an element with 0 variance! \n") ;
//   }
//   // We start the search for blocks in the lower right corner.
//   // We only do a search in columns because the Q matrix is symmetrical. Hence a search on rows would always yield
//   // the same result.
//   do {
//     nonZeros = find(vec(symmSparseMatrix(index, index, size(symmSparseMatrix.n_rows - index , 1)))) ;
//     if ((nonZeros.size() == 1)) { // We have entered a new block, this element has to be on the diagonal.
//       blockIndices.push_back(index+1) ;
//       lowerBoundary = index ;
//     }
//     index = index - 1 ;
//   } while (((max(nonZeros) + index + 1) <= lowerBoundary) & (index >= 0)) ;
//   int lastIndex = index + 2;
//   if (index < 0) {
//     lastIndex = 0 ;
//   }
//   blockIndices.push_back(lastIndex) ; // The last block will not be added in the loop, hence this step.
//   std::sort(blockIndices.begin(), blockIndices.end()) ;
//   return conv_to<uvec>::from(blockIndices) ;
// }

// void AugTree::invFromDecomposition(const sp_mat & A, const sp_mat & B, const sp_mat & D,
//                                   sp_mat * Dinv, const uvec & blockDiagAindices) {
//   cout << "Entering invFromDecomposition... \n" ;
//   sp_mat Ainv ;
//   // Typecasting those sparse matrices will be required for inversions.
//   // mat Amat = mat(A) ;
//   // mat Dmat = mat(D) ;
//   // mat Bmat = mat(B) ;
//   cout << "Inverting A... \n" ;
//   if (blockDiagAindices.size() > 1) {
//     Ainv = invertSymmBlockDiag(A, blockDiagAindices) ;
//   } else {
//     Ainv = sp_mat(inv_sympd(mat(A))) ;
//   }
//   cout << "Generating components... \n" ;
//   printf("A size: %i %i \n B size %i %i \n", A.n_rows, A.n_cols, B.n_rows, B.n_cols) ;
//   // Ainverse(0,0, size(161, 10)).print("Inverse of A:") ;
//   sp_mat M11 = sp_mat(inv(mat(A - B * (*Dinv) * trans(B)))) ;
//   cout << "Obtained M11... \n" ;
//   printf("M11 size: %i %i \n", M11.n_rows, M11.n_cols) ;
//   sp_mat M12 = -M11 * B * (*Dinv) ;
//   cout << "Obtained M12... \n" ;
//   printf("M12 size: %i %i \n", M12.n_rows, M12.n_cols) ;
//   printf("trans(B) size: %i %i \n", B.n_cols, B.n_rows) ;
//   sp_fmat identity = sp_fmat((*Dinv).n_rows, (*Dinv).n_cols) ;
//   identity.diag().ones() ;
//   sp_mat difference((*Dinv).n_rows, (*Dinv).n_cols) ;
//   sp_mat difference1 = (*Dinv) * trans(B) ;
//   difference1(0,0,size(difference1.n_rows, 1)).print("First column 1:") ;
//   sp_mat difference2 = difference1 * M12 ;
//   difference2(0,0,size(difference2.n_rows, 1)).print("First column 2:") ;
//   throw Rcpp::exception("Stop here for now...\n") ;
//   cout << "Assigned identity... \n" ;
//   cout << "Number of zeros: " << approx_equal(difference, sp_mat(difference.n_cols, difference.n_cols), "absdiff", 1e-6) << "\n";
//
//   cout << "Found difference. \n" ;
//   printf("difference size: %i %i \n", difference.n_rows, difference.n_cols) ;
//
//   // sp_mat multiplier = identity - difference ;
//   // cout << "Got multiplier... \n" ;
//   // sp_mat M22 = (*Dinv) * multiplier ;
//   // cout << "Obtained M22... \n" ;
//   // cout << "Creating the inverse matrix... \n" ;
//   // *Dinv = join_rows(join_cols(M11, trans(M12)), join_cols(M12, M22)) ; // Will the memory be freed once the function is done running?
//   // cout << "Leaving invFromDecomposition... \n" ;
// }

// std::vector<uvec> AugTree::splitBlocksDiagMatrix(std::vector<uvec> splitBlockList, uint blockSizeLimit) {
//
//   for (auto & i : splitBlockList) {
//     uint splitPos = 0 ;
//     uvec shiftedIndices = i - i.at(0) ;
//     uvec cumulativeSizes ;
//     uint cutPosition ;
//     if (i.tail(1)(0) > blockSizeLimit) {
//       cumulativeSizes = cumsum(shiftedIndices) ;
//       uvec bigPositions = find(cumulativeSizes > blockSizeLimit) ;
//       cutPosition = bigPositions.at(0) - 1 ;
//
//     }
//   }
// }

double AugTree::ComputeJointPsiMarginalPropConstant(const vec & MRAhyperStart,  const double fixedEffSDstart, const double errorSDstart, const double hyperAlpha, const double hyperBeta) {
  vec jointFunPeak = optimJointHyperMarg(MRAhyperStart, errorSDstart,fixedEffSDstart, hyperAlpha, hyperBeta) ;
  jointFunPeak.print("Optimised values:") ;
  return 0;
}

double my_f (const gsl_vector *v, void *params)  {

    gridPair *p = (gridPair *)params;
    vec MRAhyperParas(2, fill::zeros) ;
    MRAhyperParas(0)  = gsl_vector_get(v, 0) ;
    MRAhyperParas(1) = gsl_vector_get(v, 1) ;
    // The '-' is because we want to maximise. Finding the position of the minimum of "-function" is equivalent to
    // finding the position of the maximum of "function".
    return -p->grid->ComputeJointPsiMarginal(MRAhyperParas, gsl_vector_get(v, 2),
                                   gsl_vector_get(v, 3), p->vector->at(0), p->vector->at(1)) ;
  }

vec AugTree::optimJointHyperMarg(const vec & MRAhyperStart, const double fixedEffSDstart,
                  const double errorSDstart,  const double hyperAlpha, const double hyperBeta) {
  vec parasToOptim = MRAhyperStart ;
  parasToOptim.resize(MRAhyperStart.size() + 2) ;
  parasToOptim(MRAhyperStart.size()) = fixedEffSDstart ;
  parasToOptim(MRAhyperStart.size()+1) = errorSDstart ;

  uint n = parasToOptim.n_rows ;
  gridPair * par ;
  par->vector->set_size(2) ;
  par->vector->at(0) = hyperAlpha;
  par->vector->at(1) = hyperBeta;
  par->grid = this;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status;
  double size;

  /* Starting point */
  x = gsl_vector_alloc(n);

  for (int i = 0 ; i < n ; i++) {
    gsl_vector_set(x, i, parasToOptim.at(i));
  }

  /* Set initial step sizes to 0.5 */
  ss = gsl_vector_alloc(n);
  gsl_vector_set_all(ss, 0.5);

  /* Initialize method and iterate */
  minex_func.n = n ;
  minex_func.f = my_f ;
  minex_func.params = (void *)par;

  s = gsl_multimin_fminimizer_alloc(T, n);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  do
  {
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if (status)
      break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, 1e-2);

    if (status == GSL_SUCCESS)
    {
      printf ("converged to minimum at\n");
    }
  }
  while (status == GSL_CONTINUE && iter < 500);

  vec optimalParas(x->size, fill::zeros) ;

  for (uint i = 0 ; i < x->size ; i++) {
    optimalParas(i) = gsl_vector_get(x, i) ;
  }

  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (s);

  return optimalParas ;
}
