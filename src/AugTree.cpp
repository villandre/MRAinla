#include <math.h>
#include <omp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"

using namespace Rcpp ;
using namespace arma ;
using namespace MRAinla ;
using namespace Eigen ;

struct gridPair{
  AugTree * grid ;
  vec vector ;
  gridPair() { } ;
  gridPair(AugTree * gridArg, vec vectorArg) : grid(gridArg), vector(vectorArg) { } ;
};

AugTree::AugTree(uint & M, fvec & lonRange, fvec & latRange, fvec & timeRange, vec & observations, fmat & obsSp, fvec & obsTime, uint & minObsForTimeSplit, unsigned long int & seed, fmat & covariates)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime, covariates) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  m_MRAcovParas.resize(1) ;
  m_MRAcovParas.fill(-1.5e6) ;
  // m_newHyperParas = true ;

  gsl_rng_set(m_randomNumGenerator, seed) ;
  m_fixedEffParameters.resize(m_dataset.covariateValues.n_cols + 1) ; // An extra 1 is added to take into account the intercept.
  m_fixedEffParameters.randu() ;

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
  fmat spMediansMat = median(m_dataset.spatialCoords.rows(obsForMedian), 0) ;
  Row<float> spMedians = spMediansMat.row(0) ;

  fvec sortedTimes = sort(m_dataset.timeCoords.elem(obsForMedian)) ;
  float timeMedian = sortedTimes.at(std::ceil((obsForMedian.size()-1)/2));

  std::vector<dimensions> dimensionsForChildren ;
  fvec lonRange(2, fill::zeros), latRange(2, fill::zeros) ;

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
    fvec timeRange(2, fill::zeros) ;
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

void AugTree::ComputeMRAlogLik(const bool WmatsAvailable)
{
  if (!WmatsAvailable) {
    computeWmats() ;
  }
  deriveAtildeMatrices() ;

  computeOmegas() ;
  computeU() ;
  computeD() ;

  // The log-likelihood is a function of d and u computed at the root node, being at the
  // head of m_vertexVector, since it's the first one we ever created.

  double tempLogLik = (m_vertexVector.at(0)->GetD() + m_vertexVector.at(0)->GetU()) ;

  tempLogLik = -tempLogLik/2 ;
  // We set m_newHyperParas=false to indicate that covariance computations with the current hyperparameters
  // have been perfomed, and need not be performed again, unless we change the hyperparameters.
  m_recomputeMRAlogLik = false ;
  m_MRAlogLik = tempLogLik ;
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

void AugTree::computeOmegas() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveOmega(m_MRAvalues) ;
    }
  }
}

void AugTree::computeU() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = getLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveU(m_MRAvalues) ;
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
  fvec intercept = ones<fvec>(m_dataset.covariateValues.n_rows) ;
  fmat incrementedCovar = m_dataset.covariateValues ;
  incrementedCovar.insert_cols(0, intercept) ;
  vec meanVec = vec(m_dataset.covariateValues.n_rows) ;
  meanVec = incrementedCovar * m_fixedEffParameters ;
  m_dataset.responseValues -= meanVec ;
}

arma::sp_mat AugTree::CombineKandFEmatrices() {
  std::vector<mat *> KandFEmatrixList = getKmatricesInversePointers() ;
  mat FEmatrix = pow(m_fixedEffSD, 2) * eye<mat>(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  KandFEmatrixList.insert(KandFEmatrixList.begin(), &FEmatrix) ;

  // This is because we want the Q matrix to have a block-diagonal component in the lower-right corner,
  // which prompted us to make H* = [X,F] rather than H* = [F, X], like in Eidsvik.

  sp_mat Kmatrices = createSparseMatrix(KandFEmatrixList) ;

  return Kmatrices ;
}

arma::sp_mat AugTree::createHstar() {
  sp_mat Fmatrix = createFmatrix() ;

  fvec intercept = ones<fvec>(m_dataset.covariateValues.n_rows) ;
  fmat covarCopy = m_dataset.covariateValues ;
  covarCopy.insert_cols(0, intercept) ;
  sp_fmat incrementedCovar = conv_to<sp_fmat>::from(covarCopy) ;
  // We revert the order of Hstar...
  sp_mat Hstar = join_rows(conv_to<sp_mat>::from(incrementedCovar), Fmatrix) ;
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

void AugTree::ComputeLogPriors() {
  vec hyperPriorVec = m_MRAcovParas ;
  hyperPriorVec.resize(hyperPriorVec.size() + 2) ;
  hyperPriorVec(m_MRAcovParas.size()) = m_fixedEffSD ;
  hyperPriorVec(m_MRAcovParas.size() + 1) = m_errorSD ;

  std::vector<IGhyperParas> incrementedIG = m_MRAcovParasIGalphaBeta ;
  incrementedIG.push_back(m_fixedEffIGalphaBeta) ;
  incrementedIG.push_back(m_errorIGalphaBeta) ;

  hyperPriorVec.print("Hyperparameter values:") ;

  double logPrior = 0 ;

  for (uint i = 0 ; i < hyperPriorVec.size() ; i++) {
    // logPrior += hyperAlpha * log(hyperBeta) - gsl_sf_lngamma(hyperAlpha) -
    //   (hyperAlpha + 1) * log(i) - hyperBeta/i ;
    // logPrior += -((incrementedIG.at(i).m_alpha + 1) * log(hyperPriorVec.at(i)) + incrementedIG.at(i).m_beta/hyperPriorVec.at(i)) ; // Since hyperAlpha and hyperBeta do not vary, we can ignore them for optimisation purposes.
    logPrior += (incrementedIG.at(i).m_alpha - 1) * log(hyperPriorVec.at(i)) - incrementedIG.at(i).m_beta * hyperPriorVec.at(i) ;
  }

  m_logPrior = logPrior ;
}

void AugTree::ComputeLogJointCondTheta() {
  if (m_recomputeMRAlogLik) {
    ComputeMRAlogLik(true) ; // Getting the eta values requires computing the K matrices, hence KmatricesAvailable = true.
  }
  double fixedEffMean = 0 ;
  double fixedEffLogLik = 0 ;
  for (auto & i : m_fixedEffParameters) {
    fixedEffLogLik += log(normpdf(i, fixedEffMean, m_fixedEffSD)) ;
  }
  m_logCondDist = m_MRAlogLik + fixedEffLogLik ;
}


void AugTree::ComputeLogFullConditional() {

  uint n = m_dataset.responseValues.size() ;

  // sp_mat vStar = conv_to<sp_mat>::from(join_cols(m_fixedEffParameters, m_MRAetaValues)) ;

  sp_mat SigmaFEandEta = CombineKandFEmatrices() ;

  sp_mat SigmaFEandEtaInv = invertSymmBlockDiag(SigmaFEandEta, extractBlockIndices(SigmaFEandEta)) ;

  // sp_mat Hstar = createHstar() ;
  sp_mat Fmatrix = createFmatrix() ;

  sp_fmat incrementedCovar = conv_to<sp_fmat>::from(join_rows(ones<fvec>(n), m_dataset.covariateValues)) ;
  // We revert the order of Hstar...
  sp_mat Hstar = join_rows(conv_to<sp_mat>::from(incrementedCovar), Fmatrix) ;

  sp_mat TepsilonInverse = 1/pow(m_errorSD, 2) * eye<sp_mat>(n, n) ;

  sp_mat Qmat = SigmaFEandEtaInv + trans(Hstar) * TepsilonInverse * Hstar ;

  vec bVec = trans(trans(m_dataset.responseValues) * TepsilonInverse * Hstar) ;

  vec updatedMean = ComputeFullConditionalMean(bVec, Qmat) ;
  m_Vstar = updatedMean ;
  m_MRAvalues = Fmatrix * updatedMean.tail(m_numKnots) ;
  vec fixedEffMeans = updatedMean.head(m_fixedEffParameters.size()) ;
  SetFixedEffParameters(fixedEffMeans) ;

  double logDetQmat = logDeterminantQmat(Qmat) ;

  // vec recenteredVstar = m_Vstar - updatedMean ;
  // mat distExponential = trans(recenteredVstar) * Qmat * recenteredVstar ;
  // double exponential = -0.5 * distExponential.at(0, 0) ;

  // printf("Log-determinant: %.4e \n Exponent contribution: %.4e \n", detQmat, exponential) ;
  // m_logFullCond = 0.5 * detQmat + exponential ;
  m_logFullCond = 0.5 * logDetQmat ; // Since we arbitrarily evaluate always at the full-conditional mean, the exponential part of the distribution reduces to 0.
}

void AugTree::ComputeGlobalLogLik() {
  fmat incrementedCovar = m_dataset.covariateValues ;
  incrementedCovar.insert_cols(0, 1) ;
  incrementedCovar.col(0).fill(1) ;
  vec meanVec(m_dataset.responseValues.size(), fill::zeros) ;
  meanVec = incrementedCovar * m_fixedEffParameters + m_MRAvalues ;

  vec sdVec(m_dataset.responseValues.size(), fill::zeros) ;

  sdVec.fill(m_errorSD) ;
  vec logDensVec(m_dataset.responseValues.size(), fill::zeros) ;
  // logDensVec = log(normpdf(m_dataset.responseValues, meanVec, sdVec)) ;
  double logDensity = logNormPDF(m_dataset.responseValues, meanVec, sdVec) ;
  m_recomputeGlobalLogLik = false ;
  m_globalLogLik = logDensity ;
}

double AugTree::ComputeLogJointPsiMarginal() {
  // In theory, this function should not depend on the theta values...
  // We can therefore arbitrarily set them all to 0.
  cout << "Computing log-prior... \n" ;
  ComputeLogPriors() ;
  if (m_recomputeMRAlogLik) {
    computeWmats() ; // This will produce the K matrices required. NOTE: ADD CHECK THAT ENSURES THAT THE MRA LIK. IS ONLY RE-COMPUTED WHEN THE MRA COV. PARAMETERS CHANGE.
  }
  cout << "Computing log-full conditional... \n" ;
  ComputeLogFullConditional() ;
  // if (m_MRAetaValues.size() == 0) {
  //   cout << "Setting etas... \n" ;
  //   vec correlatedEtas(m_numKnots, fill::zeros) ;
  //   uint index = 0 ;
  //   arma_rng::set_seed(2) ;
  //   for (auto & i: m_vertexVector) {
  //     mat cholDecomp = chol(i->GetKmatrix()) ;
  //     vec rnormValues(cholDecomp.n_rows) ;
  //     rnormValues.randn() ;
  //     correlatedEtas.subvec(index, index + rnormValues.size()-1) =  cholDecomp * rnormValues ;
  //     index += rnormValues.size() ;
  //   }
  //   m_MRAetaValues = correlatedEtas ;
  //   sp_mat Fmatrix = createFmatrix() ;
  //   m_MRArandomValues = Fmatrix * correlatedEtas; // Pretty sure the etas match the order of knots pre-supposed by the F matrix, but would be better to make sure.
  // }
  ComputeLogJointCondTheta() ;
  uint n = m_dataset.responseValues.size() ;
  vec MRAvaluesAtObservations(m_dataset.timeCoords.n_rows) ;
  if (m_recomputeGlobalLogLik) {
    cout << "Computing log-likelihood... \n" ;
    ComputeGlobalLogLik() ;
  }
  cout << "Wrapping up... \n" ;
  printf("Total log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
         m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  return ( m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond) ;
}

// This inversion is based on recursive partitioning of the Q matrix. It is based on the observation that it is
// possible to form block-diagonal matrices on the diagonal which can be easily inverted.
// The challenge in inverting Qmat was the very heavy memory burden.
// This function involves much smaller matrices, which will make the operations easier to handle.
double AugTree::logDeterminantQmat(const sp_mat & Qmat) {
  uvec DmatrixBlockIndices = extractBlockIndicesFromLowerRight(Qmat) ;

  uint numRowsD = DmatrixBlockIndices.tail(1)(0) - DmatrixBlockIndices(0) ;
  uint shift = DmatrixBlockIndices.at(0) ;

  uvec shiftedBlockIndices = DmatrixBlockIndices - shift ;
  sp_mat Dinv = invertSymmBlockDiag(Qmat(DmatrixBlockIndices.at(0),
                                         DmatrixBlockIndices.at(0),
                                         size(numRowsD, numRowsD)),
                                         shiftedBlockIndices) ;
  double logDeterminantD = 0 ;
  for (uint i = 0 ; i < (DmatrixBlockIndices.size() - 1) ; i++) {
    double value = 0 ;
    double sign = 0 ;
    uint matSize = DmatrixBlockIndices.at(i+1) - DmatrixBlockIndices.at(i) ;
    log_det(value, sign, mat(Qmat(DmatrixBlockIndices.at(i), DmatrixBlockIndices.at(i), size(matSize, matSize)))) ;
    logDeterminantD += sign*value ;
  }
  double logDeterminantComposite, sign1 ;
  uint AmatrixSize = DmatrixBlockIndices.at(0) ;
  sp_mat compositeMat = Qmat(0, 0, size(AmatrixSize, AmatrixSize)) -
    Qmat(0, AmatrixSize, size(AmatrixSize, Qmat.n_cols - AmatrixSize)) * Dinv * Qmat(AmatrixSize, 0, size(Qmat.n_rows - AmatrixSize, AmatrixSize)) ;

  log_det(logDeterminantComposite, sign1, mat(compositeMat)) ;

  return logDeterminantD + sign1*logDeterminantComposite ;
}

uvec AugTree::extractBlockIndicesFromLowerRight(const arma::sp_mat & symmSparseMatrix) {
  std::vector<uint> blockIndices ;
  blockIndices.push_back(symmSparseMatrix.n_rows) ; // We start by adding a bound on the right, although other points in the vector correspond to upper-left corners.
  int index = symmSparseMatrix.n_rows - 2 ;
  uint lowerBoundary = symmSparseMatrix.n_rows - 1 ;
  uvec nonZeros ;
  if (any(vec(symmSparseMatrix.diag()) == 0)) {
    throw Rcpp::exception("Full conditional has an element with 0 variance! \n") ;
  }
  // We start the search for blocks in the lower right corner.
  // We only do a search in columns because the Q matrix is symmetrical. Hence a search on rows would always yield
  // the same result.
  do {
    nonZeros = find(vec(symmSparseMatrix(index, index, size(symmSparseMatrix.n_rows - index , 1)))) ;
    if ((nonZeros.size() == 1)) { // We have entered a new block, this element has to be on the diagonal.
      blockIndices.push_back(index+1) ;
      lowerBoundary = index ;
    }
    index = index - 1 ;
  } while (((max(nonZeros) + index + 1) <= lowerBoundary) & (index >= 0)) ;
  int lastIndex = index + 2;
  if (index < 0) {
    lastIndex = 0 ;
  }
  blockIndices.push_back(lastIndex) ; // The last block will not be added in the loop, hence this step.
  std::sort(blockIndices.begin(), blockIndices.end()) ;
  return conv_to<uvec>::from(blockIndices) ;
}

double AugTree::logDeterminantFullConditional(const sp_mat & SigmaMat) {
  uint n = m_dataset.timeCoords.size() ;
  sp_mat compositeAmatrix = join_rows(eye<sp_mat>(n,n),
    join_rows(conv_to<sp_mat>::from(vec(n, fill::ones)),
      conv_to<sp_mat>::from(conv_to<mat>::from(m_dataset.covariateValues)))) ;

  sp_mat Bmatrix = 1/pow(m_errorSD,2) * trans(compositeAmatrix) * eye<sp_mat>(n,n) * compositeAmatrix ;
  sp_mat Cinverse = SigmaMat ;
  sp_mat Bk(Bmatrix.n_rows , Bmatrix.n_cols) ;
  printf("Cinverse size: %i %i \n", SigmaMat.n_rows, SigmaMat.n_cols) ;
  cout << "Entering loop... \n" ;
  for (uint i = 0 ; i < 4 ; i++) {
    printf("Processing iteration %i. \n", i) ;
    uint j = Bmatrix.n_cols - i - 1;
    Bk.col(j) = Bmatrix.col(j) ;
    // double gk = 1/(1+ trace(Cinverse * trans(BkT))) ;
    double gk = 1;

    Cinverse = Cinverse - gk * Cinverse * Bk * Cinverse ;
    Bk.zeros() ;
  } ;
  throw Rcpp::exception("Stop now!!! \n") ;
  return logDeterminantQmat(Cinverse) ;
}

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

// double AugTree::ComputeJointPsiMarginalPropConstant(const vec & MRAhyperStart,
//                                 const double fixedEffSDstart, const double errorSDstart,
//                                 const double hyperAlpha, const double hyperBeta) {
//   cout << "Entered ComputeJointPsiMarginalPropConstant... \n" ;
//   vec jointFunPeak = optimJointHyperMarg(MRAhyperStart, errorSDstart, fixedEffSDstart, hyperAlpha, hyperBeta) ;
//   jointFunPeak.print("Optimised values:") ;
//   return 0;
// }

// This is for running Nelder-Mead directly in C++. It's not working yet. If optim in R doesn't deliver,
// Uncomment and fix the following.

// double my_f (const gsl_vector *v, void *params)  {
//
//     gridPair *p = (gridPair *)params;
//     vec MRAhyperParas(2, fill::zeros) ;
//     MRAhyperParas(0)  = gsl_vector_get(v, 0) ;
//     MRAhyperParas(1) = gsl_vector_get(v, 1) ;
//     // The '-' is because we want to maximise. Finding the position of the minimum of "-function" is equivalent to
//     // finding the position of the maximum of "function".
//     return -p->grid->ComputeJointPsiMarginal() ;
//   }
//
// vec AugTree::optimJointHyperMarg(const vec & MRAhyperStart, const double fixedEffSDstart,
//                   const double errorSDstart,  const double hyperAlpha, const double hyperBeta) {
//   cout << "Entered optimJointHyperMarg... \n" ;
//   vec parasToOptim = MRAhyperStart ;
//   parasToOptim.resize(MRAhyperStart.size() + 2) ;
//   parasToOptim(MRAhyperStart.size()) = fixedEffSDstart ;
//   parasToOptim(MRAhyperStart.size()+1) = errorSDstart ;
//   parasToOptim = log(parasToOptim) ;
//
//   uint n = parasToOptim.n_rows ;
//   gridPair par ;
//
//   par.vector.set_size(2) ;
//   par.vector.at(0) = hyperAlpha;
//   par.vector.at(1) = hyperBeta;
//   par.grid = this;
//   gridPair * parPoint = &par ;
//   cout << "Initializing mimimizer... \n" ;
//   const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
//   gsl_multimin_fminimizer *s = NULL;
//   gsl_vector *ss, *x;
//   gsl_multimin_function minex_func;
//
//   size_t iter = 0;
//   int status;
//   double size;
//
//   /* Starting point */
//   x = gsl_vector_alloc(n);
//
//   for (int i = 0 ; i < n ; i++) {
//     gsl_vector_set(x, i, parasToOptim.at(i));
//   }
//
//   /* Set initial step sizes to 0.5 */
//   ss = gsl_vector_alloc(n);
//   gsl_vector_set_all(ss, 0.5);
//
//   /* Initialize method and iterate */
//   minex_func.n = n ;
//   minex_func.f = my_f ;
//   minex_func.params = (void *)parPoint ;
//
//   s = gsl_multimin_fminimizer_alloc(T, n);
//   gsl_multimin_fminimizer_set(s, &minex_func, x, ss);
//   cout << "Entering minimisation loop... \n" ;
//   // do
//   // {
//   //   iter++;
//   //   status = gsl_multimin_fminimizer_iterate(s);
//   //
//   //   if (status)
//   //     break;
//   //
//   //   size = gsl_multimin_fminimizer_size(s);
//   //   status = gsl_multimin_test_size(size, 1e-2);
//   //
//   //   if (status == GSL_SUCCESS)
//   //   {
//   //     printf ("converged to minimum at\n");
//   //   }
//   // }
//   // while (status == GSL_CONTINUE && iter < 500);
//   cout << "Minimisation complete! \n" ;
//   vec optimalParas(x->size, fill::zeros) ;
//
//   for (uint i = 0 ; i < x->size ; i++) {
//     optimalParas(i) = gsl_vector_get(x, i) ;
//   }
//
//   gsl_vector_free(x);
//   gsl_vector_free(ss);
//   gsl_multimin_fminimizer_free (s);
//   cout << "Returning values... \n" ;
//   return optimalParas ;
// }

arma::vec AugTree::ComputeFullConditionalMean(const arma::vec & bVec, const arma::sp_mat & Qmat) {

  uvec blockIndices = extractBlockIndicesFromLowerRight(Qmat) ;

  uint DblockLeftIndex = blockIndices.at(0) ;
  uint sizeD = Qmat.n_cols - DblockLeftIndex + 1 ;
  uint sizeA = Qmat.n_cols - sizeD ;
  vec b1 = bVec.subvec(0, sizeA - 1) ;
  vec b2 = bVec.subvec(sizeA, bVec.size() - 1) ;
  sp_mat Bmatrix = Qmat(0, sizeA, size(sizeA, sizeD)) ;
  sp_mat Amatrix = Qmat(0, 0, size(sizeA, sizeA)) ;

  sp_mat Dinverted = invertSymmBlockDiag(Qmat(sizeA, sizeA, size(sizeD, sizeD)), blockIndices) ;

  mat secondTermInside = conv_to<mat>::from(Bmatrix * Dinverted * trans(Bmatrix)) ;

  mat compositeInverted = inv(Amatrix - secondTermInside) ;

  sp_mat firstElementSecondTerm = Dinverted * conv_to<sp_mat>::from(b2) ;
  firstElementSecondTerm = Bmatrix * firstElementSecondTerm ;
  firstElementSecondTerm = compositeInverted * firstElementSecondTerm ;
  vec firstElementFirstTerm = compositeInverted * b1 ;
  vec firstElement =  firstElementFirstTerm - firstElementSecondTerm ;

  vec secondElementFirstTerm = trans(Bmatrix) * firstElementFirstTerm ;
  secondElementFirstTerm = -Dinverted * secondElementFirstTerm ;
  sp_mat secondElementSecondTerm = trans(Bmatrix) * firstElementSecondTerm ;
  secondElementSecondTerm = Dinverted * secondElementSecondTerm ;
  secondElementSecondTerm = Dinverted * b2 + secondElementSecondTerm ;
  vec secondElement = secondElementFirstTerm + secondElementSecondTerm ;

  vec meanVec = join_cols<vec>(firstElement, secondElement) ;

  return meanVec ;
}

