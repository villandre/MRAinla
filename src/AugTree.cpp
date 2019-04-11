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

struct gridPair{
  AugTree * grid ;
  vec vector ;
  gridPair() { } ;
  gridPair(AugTree * gridArg, vec vectorArg) : grid(gridArg), vector(vectorArg) { } ;
};

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, vec & timeRange, vec & observations, mat & obsSp, vec & obsTime, uint & minObsForTimeSplit, unsigned long int & seed, mat & covariates)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime, covariates) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  m_MRAcovParas.resize(1) ;
  m_MRAcovParas.fill(-1.5e6) ;

  gsl_rng_set(m_randomNumGenerator, seed) ;
  m_fixedEffParameters.resize(m_dataset.covariateValues.n_cols + 1) ; // An extra 1 is added to take into account the intercept.
  m_fixedEffParameters.randu() ;

  BuildTree(minObsForTimeSplit) ;

  m_numKnots = 0 ;
  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    m_vertexVector.at(i)->SetNodeId(i) ;
    m_numKnots += m_vertexVector.at(i)->GetKnotsCoor().timeCoords.size() ;
  }
  m_FullCondMean = vec(m_numKnots + m_fixedEffParameters.size(), fill::zeros) ;
  m_FullCondSDs = vec(m_numKnots + m_fixedEffParameters.size(), fill::zeros) ;
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
  Row<double> spMedians = spMediansMat.row(0) ;

  vec sortedTimes = sort(m_dataset.timeCoords.elem(obsForMedian)) ;
  double timeMedian = sortedTimes.at(std::ceil((obsForMedian.size()-1)/2));

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
    vec timeRange(2, fill::zeros) ;
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
  int incrementedDepth = parent->GetDepth()+1 ;

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
  // double scaledRespSize =  0.05 ;
  // double depthContribution = 19*node->GetDepth() ;
  // depthContribution = depthContribution/60 ;
  // double expToRound = depthContribution + scaledRespSize ; // Why is it that we have to store every term separately for expToRound to not evaluate to 0?
  // double numKnots = expToRound * node->GetObsInNode().size() ;
  //// CAREFUL: THIS IS HARD-CODED, THERE MIGHT BE A BETTER WAY /////////////
  int baseNumberOfKnots = 30 ;
  uint numKnotsToGen = baseNumberOfKnots * pow(2, node->GetDepth()) ;
  // double numKnots = sqrt((double) node->GetObsInNode().size()) ;
  // uint numKnotsToGen = uint (std::ceil(numKnots)) ;
  node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i) ;
    }
  }
}

// void AugTree::ComputeMRAlogLik(const bool WmatsAvailable)
// {
//   cout << "Entered ComputeMRAlogLik... \n" ;
//   if (!WmatsAvailable) {
//     computeWmats() ;
//   }
//   deriveAtildeMatrices() ;
//
//   computeOmegas() ;
//   computeU() ;
//   computeD() ;
//
//   // The log-likelihood is a function of d and u computed at the root node, being at the
//   // head of m_vertexVector, since it's the first one we ever created.
//
//   double tempLogLik = (m_vertexVector.at(0)->GetD() + m_vertexVector.at(0)->GetU()) ;
//
//   tempLogLik = -tempLogLik/2 ;
//   // We set m_newHyperParas=false to indicate that covariance computations with the current hyperparameters
//   // have been perfomed, and need not be performed again, unless we change the hyperparameters.
//   m_recomputeMRAlogLik = false ;
//   m_MRAlogLik = tempLogLik ;
// }

void AugTree::ComputeMRAlogLikAlt(const bool WmatsAvailable)
{
  cout << "Entered ComputeMRAlogLik... \n" ;
  if (!WmatsAvailable) {
    computeWmats() ;
  }
  double MRAlogLik = 0 ;
  int currentIndex = m_fixedEffParameters.size() ;

  for (auto & i : m_vertexVector) {
    int lastIndex = currentIndex + i->GetNumKnots() - 1 ;
    vec etas = m_Vstar.subvec(currentIndex, lastIndex) ;
    // etas.print("Eta vector:") ;

    double logDeterminantValue = 0 ;
    double sign = 0 ;
    log_det(logDeterminantValue, sign, i->GetKmatrix()) ;
    if (sign < 0) {
      throw Rcpp::exception("Error in ComputeMRAlogLikAlt! sign should be positive. \n") ;
    }

    mat exponentTerm = -0.5 * trans(etas) * i->GetKmatrixInverse() * etas ;
    MRAlogLik += (-0.5*logDeterminantValue + exponentTerm(0,0)) ;
    currentIndex = lastIndex + 1 ;
  }

  m_MRAlogLik = MRAlogLik ;
}

void AugTree::computeWmats() {

  m_vertexVector.at(0)->ComputeWmat(m_MRAcovParas, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;

  for (uint level = 1; level <= m_M; level++) {
    std::vector<TreeNode *> levelNodes = GetLevelNodes(level) ;

    // Trying openmp. We need to have a standard looping structure.
    // #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++)
    {
      (*it)->ComputeWmat(m_MRAcovParas, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;
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

void AugTree::deriveAtildeMatrices() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;

    // for (auto & i : levelNodes) {
    //   i->DeriveAtilde() ;
    // }
    // #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
      (*it)->DeriveAtilde() ;
    }
  }
}

void AugTree::computeOmegas() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveOmega(m_MRAvalues) ;
    }
  }
}

void AugTree::computeU() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;
    for (auto & i : levelNodes) {
      i->DeriveU(m_MRAvalues) ;
    }
  }
};

void AugTree::computeD() {
  for (int level = m_M; level >= 0 ; level--) {
    uint levelRecast = (uint) level ;
    std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;

    // Trying openmp. We need to have a standard looping structure.

    // #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
      (*it)->DeriveD() ;
    }
  }
}

// std::vector<GaussDistParas> AugTree::ComputeConditionalPrediction(const inputdata & predictionData) {
//   distributePredictionData(predictionData) ;
//
//   for (int level = m_M; level >= 0 ; level--) {
//     uint levelRecast = (uint) level ;
//     std::vector<TreeNode *> levelNodes = GetLevel(levelRecast) ;
//     // #pragma omp parallel for
//     for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
//     // for (auto & i : levelNodes) {
//       (*it)->ComputeParasEtaDeltaTilde(predictionData, m_dataset, m_MRAcovParas) ;
//     }
//   }
//   computeBtildeInTips() ;
//
//   std::vector<vec> predictionsFromEachSection ;
//   std::vector<TreeNode *> tipNodes = GetLevel(m_M) ;
//
//   std::vector<GaussDistParas> distParasFromEachZone ;
//   // #pragma omp parallel for
//   for (std::vector<TreeNode *>::iterator it = tipNodes.begin(); it < tipNodes.end(); it++) {
//   // for (auto & i : tipNodes) {
//     distParasFromEachZone.push_back((*it)->CombineEtaDelta(predictionData, m_fixedEffParameters)) ;
//   }
//
//   return distParasFromEachZone ;
// }


//
// void AugTree::computeBtildeInTips() {
//   std::vector<std::vector<TreeNode *>> nodesAtLevels(m_M+1) ;
//   for (uint i = 0 ; i <= m_M ; i++) {
//     nodesAtLevels.at(i) = GetLevel(i) ;
//   }
//   // The two following loops cannot be joined.
//   for (uint level = 1; level <= m_M ; level++) {
//     for (auto & i : nodesAtLevels.at(level)) {
//       i->initiateBknots(m_MRAcovParas) ;
//     }
//   }
//   for (uint level = 1; level < m_M ; level++) {
//     for (auto & i : nodesAtLevels.at(level+1)) {
//       i->completeBknots(m_MRAcovParas, level) ;
//     }
//   }
//   for (auto & i : nodesAtLevels.at(m_M)) {
//     i->computeBpred(m_predictLocations, m_MRAcovParas) ;
//     i->deriveBtilde(m_predictLocations) ;
//   }
// }

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

arma::sp_mat AugTree::CombineFEinvAndKinvMatrices() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  mat FEinvMatrix = pow(m_fixedEffSD, -2) * eye<mat>(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;

  // This is because we want the Q matrix to have a block-diagonal component in the lower-right corner,
  // which prompted us to make H* = [X,F] rather than H* = [F, X], like in Eidsvik.

  sp_mat FEinvAndKinvMatrices = createBlockMatrix(FEinvAndKinvMatrixList) ;

  return FEinvAndKinvMatrices ;
}

// arma::sp_mat AugTree::createHstar() {
//   sp_mat Fmatrix = createFmatrixAlt() ;
//
//   vec intercept = ones<vec>(m_dataset.covariateValues.n_rows) ;
//   mat covarCopy = m_dataset.covariateValues ;
//   covarCopy.insert_cols(0, intercept) ;
//   sp_mat incrementedCovar = conv_to<sp_mat>::from(covarCopy) ;
//   sp_mat Hstar = join_rows(conv_to<sp_mat>::from(incrementedCovar), Fmatrix) ;
//   return Hstar ;
// }

arma::sp_mat AugTree::createFmatrixAlt(const bool predictionMode) {

  int numObs = m_dataset.spatialCoords.n_rows ;
  if (predictionMode) {
    numObs = m_predictData.spatialCoords.n_rows ;
  }
  sp_mat Fmat(numObs, m_numKnots) ;
  std::vector<uvec> FmatNodeOrder(m_M) ;
  uvec FmatObsOrder(numObs, fill::zeros) ;
  uint rowIndex = 0 ;
  uvec colIndicesPerNodeId(m_vertexVector.size(), fill::zeros) ;

  uint colIndex = 0 ;
  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    colIndicesPerNodeId.at(i) = colIndex ;
    colIndex += m_vertexVector.at(i)->GetNumKnots() ;
  }

  std::vector<TreeNode *> tipNodes = GetLevelNodes(m_M) ;
  int numTips = tipNodes.size() ;

  std::vector<uvec> ancestorIdsVec(numTips) ;

  for (uint i = 0 ; i < tipNodes.size(); i++) {
    uvec idVec = tipNodes.at(i)->GetAncestorIds() ; // Last element is tip node.
    ancestorIdsVec.at(i) = idVec ;
  }

  std::sort(tipNodes.begin(), tipNodes.end(), [] (TreeNode * first, TreeNode * second) {
    uvec firstAncestorIds = first->GetAncestorIds() ;
    uvec secondAncestorIds = second->GetAncestorIds() ;
    bool firstSmallerThanSecond = false ;
    for (uint i = 1 ; i < firstAncestorIds.size() ; i++) { // The ultimate ancestor is always node 0, hence the loop starting at 1.
      bool test = (firstAncestorIds.at(i) == secondAncestorIds.at(i)) ;
      if (!test) {
        firstSmallerThanSecond = (firstAncestorIds.at(i) < secondAncestorIds.at(i)) ;
        break ;
      }
    }
    return firstSmallerThanSecond ;
  }) ; // This is supposed to reorder the tip nodes in such a way that the F matrix has contiguous blocks.

  for (auto & nodeToProcess : tipNodes) {
    uvec observationIndices = nodeToProcess->GetObsInNode() ;

    if (predictionMode) {
      observationIndices = nodeToProcess->GetPredIndices() ;
    }
    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.subvec(rowIndex, size(observationIndices)) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      for (uint resolution = 0; resolution <= m_M; resolution++) {
        uint nodeId = brickList.at(resolution)->GetNodeId() ;
        mat matToAdd ;
        if (predictionMode) {
          matToAdd = nodeToProcess->GetUpred(resolution) ;
        } else {
          matToAdd = nodeToProcess->GetB(resolution) ;
        }

        Fmat(rowIndex, colIndicesPerNodeId.at(nodeId), size(matToAdd)) = matToAdd ;
      }
      // rowIndex += nodeToProcess->GetB(0).n_rows ; // The B(0), B(1), ..., B(M) all have the same number of rows.
      rowIndex += observationIndices.size() ; // The B matrices should have as may rows as observations in the node...
    }
  }

  m_obsOrderForFmat = FmatObsOrder ;
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
    logPrior += (incrementedIG.at(i).m_alpha - 1) * log(hyperPriorVec.at(i)) - incrementedIG.at(i).m_beta * hyperPriorVec.at(i) ;
  }

  m_logPrior = logPrior ;
}

// void AugTree::ComputeLogJointCondTheta() {
//   // if (m_recomputeMRAlogLik) {
//     ComputeMRAlogLikAlt(true) ; // Getting the eta values requires computing the K matrices, hence KmatricesAvailable = true.
//   // }
//
//   double fixedEffLogLik = 0 ;
//
//   for (uint i = 0; i < m_fixedEffParameters.size(); i++) {
//     fixedEffLogLik += std::log(normpdf(m_fixedEffParameters.at(i), m_FEmu.at(i), m_fixedEffSD)) ;
//   }
//
//   m_logCondDist = m_MRAlogLik + fixedEffLogLik ;
// }


void AugTree::ComputeLogFCandLogCDandDataLL() {

  int n = m_dataset.responseValues.size() ;

  sp_mat SigmaFEandEtaInv = CombineFEinvAndKinvMatrices() ;

  sp_mat Fmatrix = createFmatrixAlt(false) ;

  mat transIncrementedCovar = trans(join_rows(ones<vec>(n), m_dataset.covariateValues)) ;

 //We re-shuffle the X matrix in such a way that the lines match those in the F matrix, based on
 // m_obsOrderForFmat.

  mat incrementedCovarReshuffled = trans(transIncrementedCovar.cols(m_obsOrderForFmat)) ;

  sp_mat Hstar = join_rows(conv_to<sp_mat>::from(incrementedCovarReshuffled), Fmatrix) ;

  sp_mat TepsilonInverse = std::pow(m_errorSD, -2) * eye<sp_mat>(n, n) ;

  m_FullCondPrecision = SigmaFEandEtaInv + trans(Hstar) * TepsilonInverse * Hstar ;

  vec responsesReshuffled = m_dataset.responseValues.elem(m_obsOrderForFmat) ;

  vec bVec = trans(trans(responsesReshuffled) * TepsilonInverse * Hstar) ;

  vec updatedMean = ComputeFullConditionalMean(bVec) ;

  m_Vstar = updatedMean ;

  m_MRAvalues = Fmatrix * updatedMean.tail(m_numKnots) ;

  vec fixedEffMeans = updatedMean.head(m_fixedEffParameters.size()) ;
  SetFixedEffParameters(fixedEffMeans) ;

  double logDetQmat = logDeterminantQmat() ;

  // sp_mat Kmatrices = createSparseMatrix(KmatrixList) ;

  m_logFullCond = 0.5 * logDetQmat ; // Since we arbitrarily evaluate always at the full-conditional mean, the exponential part of the distribution reduces to 0.

  // Computing p(v* | Psi)
  mat logLikTerm = -0.5 * trans(m_Vstar) * SigmaFEandEtaInv * m_Vstar ;
  double logDetSigmaKFEinv = logDetBlockMatrix(SigmaFEandEtaInv) ;
  m_logCondDist = logLikTerm(0,0) + 0.5 * logDetSigmaKFEinv ;

  // Computing p(y | v*, Psi)
  double errorLogDet = -2 * n * std::log(m_errorSD) ;
  vec recenteredY = responsesReshuffled - Hstar * m_Vstar ;
  vec globalLogLikExp = -0.5 * (trans(recenteredY) * TepsilonInverse * recenteredY) ;
  m_globalLogLik = 0.5 * errorLogDet + globalLogLikExp(0) ;
}

// void AugTree::ComputeGlobalLogLik() {
//   mat incrementedCovar = m_dataset.covariateValues ;
//   incrementedCovar.insert_cols(0, 1) ;
//   incrementedCovar.col(0).fill(1) ;
//   mat incrementedCovarTrans = trans(incrementedCovar) ;
//   incrementedCovar = trans(incrementedCovarTrans.cols(m_obsOrderForFmat)) ;
//   vec meanVec(m_dataset.responseValues.size(), fill::zeros) ;
//   meanVec = incrementedCovar * m_fixedEffParameters + m_MRAvalues ;
//
//   vec sdVec(m_dataset.responseValues.size(), fill::zeros) ;
//
//   sdVec.fill(m_errorSD) ;
//   vec logDensVec(m_dataset.responseValues.size(), fill::zeros) ;
//   // logDensVec = log(normpdf(m_dataset.responseValues, meanVec, sdVec)) ;
//   double logDensity = logNormPDF(m_dataset.responseValues.elem(m_obsOrderForFmat), meanVec, sdVec) ;
//   m_recomputeGlobalLogLik = false ;
//   m_globalLogLik = logDensity ;
// }

void AugTree::ComputeLogJointPsiMarginal() {

  ComputeLogPriors() ;

  if (m_recomputeMRAlogLik) {
    computeWmats() ; // This will produce the K matrices required. NOTE: ADD CHECK THAT ENSURES THAT THE MRA LIK. IS ONLY RE-COMPUTED WHEN THE MRA COV. PARAMETERS CHANGE.
  }

  ComputeLogFCandLogCDandDataLL() ;

  printf("Observations log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
         m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  m_logJointPsiMarginal = m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond ;
}

// This inversion is based on recursive partitioning of the Q matrix. It is based on the observation that it is
// possible to form block-diagonal matrices on the diagonal which can be easily inverted.
// The challenge in inverting Qmat was the very heavy memory burden.
// This function involves much smaller matrices, which will make the operations easier to handle.
double AugTree::logDeterminantQmat() {
  uvec DmatrixBlockIndices = extractBlockIndicesFromLowerRight(m_FullCondPrecision) ;

  int numRowsD = m_FullCondPrecision.n_rows - DmatrixBlockIndices(0) ;

  int shift = DmatrixBlockIndices.at(0) ;

  uvec shiftedBlockIndices = DmatrixBlockIndices - shift ;
  sp_mat Dinv = invertSymmBlockDiag(m_FullCondPrecision(DmatrixBlockIndices.at(0),
                                         DmatrixBlockIndices.at(0),
                                         size(numRowsD, numRowsD)),
                                         shiftedBlockIndices) ;
  double logDeterminantD = 0 ;
  for (uint i = 0 ; i < (DmatrixBlockIndices.size() - 1) ; i++) {
    double value = 0 ;
    double sign = 0 ;
    uint matSize = DmatrixBlockIndices.at(i+1) - DmatrixBlockIndices.at(i) ;
    log_det(value, sign, mat(m_FullCondPrecision(DmatrixBlockIndices.at(i), DmatrixBlockIndices.at(i), size(matSize, matSize)))) ;
    if (sign < 0) {
      throw Rcpp::exception("Error in logDeterminantQmat! sign should be positive. \n") ;
    }
    logDeterminantD += value ;
  }
  double logDeterminantComposite, sign1 ;
  uint AmatrixSize = DmatrixBlockIndices.at(0) ;
  sp_mat compositeMat = m_FullCondPrecision(0, 0, size(AmatrixSize, AmatrixSize)) -
    m_FullCondPrecision(0, AmatrixSize, size(AmatrixSize, numRowsD)) * Dinv * m_FullCondPrecision(AmatrixSize, 0, size(numRowsD, AmatrixSize)) ;

  log_det(logDeterminantComposite, sign1, mat(compositeMat)) ;
  if (sign1 < 0) {
    throw Rcpp::exception("Error in logDeterminantQmat! sign1 should be positive. \n") ;
  }

  return logDeterminantD + logDeterminantComposite ;
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

// double AugTree::logDeterminantFullConditional(const sp_mat & SigmaMat) {
//   uint n = m_dataset.timeCoords.size() ;
//   sp_mat compositeAmatrix = join_rows(eye<sp_mat>(n,n),
//     join_rows(conv_to<sp_mat>::from(vec(n, fill::ones)),
//       conv_to<sp_mat>::from(conv_to<mat>::from(m_dataset.covariateValues)))) ;
//
//   sp_mat Bmatrix = 1/pow(m_errorSD,2) * trans(compositeAmatrix) * eye<sp_mat>(n,n) * compositeAmatrix ;
//   sp_mat Cinverse = SigmaMat ;
//   sp_mat Bk(Bmatrix.n_rows , Bmatrix.n_cols) ;
//
//   for (uint i = 0 ; i < 4 ; i++) {
//
//     uint j = Bmatrix.n_cols - i - 1;
//     Bk.col(j) = Bmatrix.col(j) ;
//     // double gk = 1/(1+ trace(Cinverse * trans(BkT))) ;
//     double gk = 1;
//
//     Cinverse = Cinverse - gk * Cinverse * Bk * Cinverse ;
//     Bk.zeros() ;
//   } ;
//   throw Rcpp::exception("Stop now!!! \n") ;
//   return logDeterminantQmat(Cinverse) ;
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
//   sp_mat identity = sp_mat((*Dinv).n_rows, (*Dinv).n_cols) ;
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

arma::vec AugTree::ComputeFullConditionalMean(const arma::vec & bVec) {

  uvec blockIndices = extractBlockIndicesFromLowerRight(m_FullCondPrecision) ;

  uint DblockLeftIndex = blockIndices.at(0) ;
  uint sizeD = m_FullCondPrecision.n_cols - DblockLeftIndex ;
  uint sizeA = m_FullCondPrecision.n_cols - sizeD ;
  vec b1 = bVec.subvec(0, sizeA - 1) ;
  vec b2 = bVec.subvec(sizeA, bVec.size() - 1) ;
  sp_mat Bmatrix = m_FullCondPrecision(0, sizeA, size(sizeA, sizeD)) ;
  sp_mat Amatrix = m_FullCondPrecision(0, 0, size(sizeA, sizeA)) ;

  sp_mat Dinverted = invertSymmBlockDiag(m_FullCondPrecision(sizeA, sizeA, size(sizeD, sizeD)), blockIndices) ;

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

  if (m_recordFullConditional) {
    vec firstDiag = sqrt(compositeInverted.diag()) ;
    m_FullCondSDs.subvec(0, size(firstDiag)) = firstDiag ;

    for (uint i = 0 ; i < Dinverted.n_rows ; i++) {
      sp_mat firstExp = Bmatrix * Dinverted.col(i) ;
      mat secondExp = compositeInverted * conv_to<mat>::from(firstExp) ;
      mat thirdExp = trans(Bmatrix) * secondExp ;
      mat finalExp = Dinverted.col(i) + Dinverted * thirdExp ;
      m_FullCondSDs.at(i + compositeInverted.n_rows) = std::sqrt(finalExp(i, 0)) ;
    }

    m_FullCondMean = meanVec ;
  }

  return meanVec ;
}

sp_mat AugTree::ComputeHpred(const mat & spCoords, const vec & time, const mat & covariateMatrix) {

  SetPredictData(spCoords, time, covariateMatrix) ;
  std::vector<TreeNode *> tipNodes = GetLevelNodes(m_M) ;

  for (auto & i : tipNodes) {
    i->SetPredictLocations(m_predictData) ;
    uvec predictionsInLeaf = i->GetPredIndices() ;
    if (predictionsInLeaf.size() > 0) {

      i->computeUpred(m_MRAcovParas, m_predictData, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;
    }
  }
  sp_mat FmatrixPred = createFmatrixAlt(true) ;
  mat incrementedCovariate = join_rows(vec(covariateMatrix.n_rows, fill::ones), covariateMatrix) ;
  sp_mat HmatrixPred = join_rows(conv_to<sp_mat>::from(incrementedCovariate), FmatrixPred) ;
  // conv_to<mat>::from(HmatrixPred).save("/home/luc/Documents/HmatrixPred.info", raw_ascii) ;
  return HmatrixPred ;
}

arma::vec AugTree::ComputeEvar(const arma::sp_mat & HmatPred) {
  uvec blockIndices = extractBlockIndicesFromLowerRight(m_FullCondPrecision) ;
  double errorVar = std::pow(m_errorSD, 2) ;
  vec EvarValues(HmatPred.n_rows, fill::zeros) ;

  uint DblockLeftIndex = blockIndices.at(0) ;
  uint sizeD = m_FullCondPrecision.n_cols - DblockLeftIndex ;
  uint sizeA = m_FullCondPrecision.n_cols - sizeD ;

  sp_mat Bmatrix = m_FullCondPrecision(0, sizeA, size(sizeA, sizeD)) ;
  sp_mat Amatrix = m_FullCondPrecision(0, 0, size(sizeA, sizeA)) ;
  sp_mat Dinverted = invertSymmBlockDiag(m_FullCondPrecision(sizeA, sizeA, size(sizeD, sizeD)), blockIndices) ;
  mat secondTermInside = conv_to<mat>::from(Bmatrix * Dinverted * trans(Bmatrix)) ;
  sp_mat compositeInverted = conv_to<sp_mat>::from(inv(Amatrix - secondTermInside)) ;
  // mat compositeInverted = inv(Amatrix - secondTermInside) ;

  sp_mat bVec ;
  sp_mat b1, b2 ;

  for (uint i = 0 ; i < HmatPred.n_rows ; i++) {
    bVec = trans(HmatPred.row(i)) ;
    b1 = bVec(0, 0, size(sizeA, 1)) ;
    b2 = bVec(sizeA, 0, size(sizeD, 1)) ;

    sp_mat DinvB2 = Dinverted * b2 ;
    sp_mat BmatDinvB2 = Bmatrix * DinvB2 ;
    sp_mat firstElement = b1 - BmatDinvB2 ;
    sp_mat b1CompInv = trans(compositeInverted * b1) ;
    sp_mat firstTerm =  b1CompInv * firstElement ;

    sp_mat secondTerm = (-b1 + BmatDinvB2) ;
    secondTerm = compositeInverted * secondTerm ;
    secondTerm = trans(BmatDinvB2) * secondTerm + trans(DinvB2) * b2 ;

    sp_mat meanVec = firstTerm + secondTerm  ;
    // sp_mat firstElementSecondTerm = Dinverted * b2 ;
    // firstElementSecondTerm = Bmatrix * firstElementSecondTerm ;
    // firstElementSecondTerm = compositeInverted * firstElementSecondTerm ;
    // sp_mat firstElementFirstTerm = conv_to<sp_mat>::from(compositeInverted) * b1 ;
    // sp_mat firstElement =  firstElementFirstTerm - firstElementSecondTerm ;
    //
    // sp_mat secondElementFirstTerm = trans(Bmatrix) * firstElementFirstTerm ;
    // secondElementFirstTerm = -Dinverted * secondElementFirstTerm ;
    // sp_mat secondElementSecondTerm = trans(Bmatrix) * firstElementSecondTerm ;
    // secondElementSecondTerm = Dinverted * secondElementSecondTerm ;
    // secondElementSecondTerm = Dinverted * b2 + secondElementSecondTerm ;
    // sp_mat secondElement = secondElementFirstTerm + secondElementSecondTerm ;
    //
    // sp_mat meanVec = trans(bVec) * join_cols(firstElement, secondElement) ;
    EvarValues.at(i) = meanVec.at(0,0) + errorVar ;
  }
  return EvarValues ;
}
