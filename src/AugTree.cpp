#include <math.h>
// #include <omp.h>
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

struct MVN{
  vec mu ;
  sp_mat precision ;

  MVN() { } ;
  MVN(const vec & mu, const sp_mat & precision): mu(mu), precision(precision) {} ;

  double ComputeExpTerm(const vec & coordinate) {
    vec output = -0.5 * trans(coordinate - mu) * precision * (coordinate - mu) ;
    return output(0) ;
  }
};

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, vec & timeRange, vec & observations, mat & obsSp, vec & obsTime, mat & predCovariates, mat & predSp, vec & predTime, uint & minObsForTimeSplit, unsigned long int & seed, mat & covariates, const bool splitTime, const unsigned int numKnotsRes0, const unsigned int J)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime, covariates) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  SetPredictData(predSp, predTime, predCovariates) ;

  gsl_rng_set(m_randomNumGenerator, seed) ;
  m_fixedEffParameters.resize(m_dataset.covariateValues.n_cols + 1) ; // An extra 1 is added to take into account the intercept.
  m_fixedEffParameters.randu() ;

  BuildTree(minObsForTimeSplit, splitTime, numKnotsRes0, J) ;

  m_numKnots = 0 ;

  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    m_numKnots += m_vertexVector.at(i)->GetKnotsCoor().timeCoords.size() ;
  }
  m_FullCondMean = vec(m_numKnots + m_fixedEffParameters.size(), fill::zeros) ;
  m_FullCondSDs = vec(m_numKnots + m_fixedEffParameters.size(), fill::zeros) ;
  m_Vstar = vec(m_numKnots + m_dataset.covariateValues.n_cols + 1, fill::zeros) ;
  m_Vstar.subvec(0, size(m_FEmu)) = m_FEmu ;
  m_HmatPos = umat(0, 2) ;
  // createHmatrixPos() ;
  m_HmatPredPos = umat(0, 2) ;
  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

  for (auto & i : tipNodes) {
    i->SetPredictLocations(m_predictData) ;
  }

  createHmatrixPredPos() ;
  m_SigmaPos = umat(0, 2) ;
}

void AugTree::BuildTree(const uint & minObsForTimeSplit, const bool splitTime, const unsigned int numKnots0, const unsigned int J)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, m_dataset) ;

  m_vertexVector.push_back(topNode) ;

  createLevels(topNode, minObsForTimeSplit, splitTime) ;
  numberNodes() ;

  generateKnots(topNode, numKnots0, J) ;
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

// In this version, we first split the latitude, then the latitude, and finally time.
// We make sure that splits don't result in empty regions
void AugTree::createLevels(TreeNode * parent, const uint & numObsForTimeSplit, const bool splitTime) {
  uvec obsForMedian = parent->GetObsInNode() ;
  uvec childMembership(obsForMedian.size(), fill::zeros) ;
  std::vector<dimensions> childDimensions ;
  childDimensions.push_back(parent->GetDimensions()) ;

  if (obsForMedian.n_elem <= 1) {
    throw Rcpp::exception("Cannot split empty region or region with only one observation.\n") ;
  }

  uint newChildIndex = 1 ;

  // Handling longitude
  uvec currentChildIndices = regspace<uvec>(0, newChildIndex - 1) ;
  for (auto & j : currentChildIndices) {
    uvec elementsInChild = find(childMembership == j) ;
    vec column = m_dataset.spatialCoords.col(0) ;
    vec subColumn = column.elem(obsForMedian) ;
    double colMedian = median(subColumn.elem(elementsInChild)) ;
    vec updatedLongitude{childDimensions.at(j).longitude.at(0), colMedian} ;
    vec newChildLongitude{colMedian, childDimensions.at(j).longitude.at(1)} ;
    dimensions newDimensions = childDimensions.at(j) ;
    newDimensions.longitude = newChildLongitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(j).longitude = updatedLongitude ;
    uvec greaterElements = find(subColumn.elem(elementsInChild) > colMedian) ; // In deriveObsInNode, the checks are <=. It follows that observations on the right and upper boundaries of a zone are included.
    uvec updateIndices = elementsInChild.elem(greaterElements) ;
    childMembership.elem(updateIndices).fill(newChildIndex) ;
    newChildIndex += 1 ;
  }

  // Handling latitude
  currentChildIndices = regspace<uvec>(0, newChildIndex - 1) ;
  for (auto & j : currentChildIndices) {
    vec column = m_dataset.spatialCoords.col(1) ;
    vec subColumn = column.elem(obsForMedian) ;
    uvec elementsInChild = find(childMembership == j) ;
    double colMedian = median(subColumn.elem(elementsInChild)) ;
    vec updatedLatitude{childDimensions.at(j).latitude.at(0), colMedian} ;
    vec newChildLatitude{colMedian, childDimensions.at(j).latitude.at(1)} ;
    dimensions newDimensions = childDimensions.at(j) ;
    newDimensions.latitude = newChildLatitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(j).latitude = updatedLatitude ;
    uvec greaterElements = find(subColumn.elem(elementsInChild) > colMedian) ;
    uvec updateIndices = elementsInChild.elem(greaterElements) ;
    childMembership.elem(updateIndices).fill(newChildIndex) ;
    newChildIndex += 1 ;
  }
  if (splitTime) {
    // Handling time
    currentChildIndices = regspace<uvec>(0, newChildIndex - 1) ;
    for (auto & j : currentChildIndices) {
      vec time = m_dataset.timeCoords.elem(obsForMedian) ;
      uvec elementsInChild = find(childMembership == j) ;
      vec elementsForMedian = time.elem(elementsInChild) ;
      uvec uniqueTimeValues = find_unique(elementsForMedian) ;

      if (uniqueTimeValues.size() > 1) {
        double colMedian = median(elementsForMedian) ;
        vec updatedTime{childDimensions.at(j).time.at(0), colMedian} ;
        vec newChildTime{colMedian, childDimensions.at(j).time.at(1)} ;
        dimensions newDimensions = childDimensions.at(j) ;
        newDimensions.time = newChildTime ;
        childDimensions.push_back(newDimensions) ;
        childDimensions.at(j).time = updatedTime ;
        uvec greaterElements = find(time.elem(elementsInChild) > colMedian) ;
        uvec updateIndices = elementsInChild.elem(greaterElements) ;
        childMembership.elem(updateIndices).fill(newChildIndex) ;
        newChildIndex += 1 ;
      }
    }
  }

  int incrementedDepth = parent->GetDepth()+1 ;
  for (auto & i : childDimensions) {
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
      createLevels(i, numObsForTimeSplit, splitTime) ;
    }
  }
}

void AugTree::generateKnots(TreeNode * node, const unsigned int numKnotsRes0, const unsigned int J) {

  int numNodesAtLevel = GetLevelNodes(node->GetDepth()).size() ;
  uint numKnotsToGen = std::max(uint(std::ceil((numKnotsRes0 * pow(J, node->GetDepth()))/numNodesAtLevel)), uint(2)) ;

  node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;

  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i, numKnotsRes0, J) ;
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
  if (!WmatsAvailable) {
    computeWmats() ;
  }
  double MRAlogLik = 0 ;
  int currentIndex = m_fixedEffParameters.size() ;

  for (auto & i : m_vertexVector) {
    int lastIndex = currentIndex + i->GetNumKnots() - 1 ;
    vec etas = m_Vstar.subvec(currentIndex, lastIndex) ;

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

  m_vertexVector.at(0)->ComputeWmat(m_MRAcovParasSpace, m_MRAcovParasTime, m_spacetimeScaling, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;

  for (uint level = 1; level <= m_M; level++) {
    std::vector<TreeNode *> levelNodes = GetLevelNodes(level) ;

    // Trying openmp. We need to have a standard looping structure.
    // #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++)
    {
      (*it)->ComputeWmat(m_MRAcovParasSpace, m_MRAcovParasTime, m_spacetimeScaling, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;
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

// void AugTree::computeOmegas() {
//   for (int level = m_M; level >= 0 ; level--) {
//     uint levelRecast = (uint) level ;
//     std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;
//     for (auto & i : levelNodes) {
//       i->DeriveOmega(m_MRAvalues) ;
//     }
//   }
// }
//
// void AugTree::computeU() {
//   for (int level = m_M; level >= 0 ; level--) {
//     uint levelRecast = (uint) level ;
//     std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;
//     for (auto & i : levelNodes) {
//       i->DeriveU(m_MRAvalues) ;
//     }
//   }
// };
//
// void AugTree::computeD() {
//   for (int level = m_M; level >= 0 ; level--) {
//     uint levelRecast = (uint) level ;
//     std::vector<TreeNode *> levelNodes = GetLevelNodes(levelRecast) ;
//
//     // Trying openmp. We need to have a standard looping structure.
//
//     // #pragma omp parallel for
//     for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
//       (*it)->DeriveD() ;
//     }
//   }
// }

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

arma::sp_mat AugTree::CreateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  mat FEinvMatrix = pow(m_fixedEffSD, -2) * eye<mat>(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;

  // This is because we want the Q matrix to have a block-diagonal component in the lower-right corner,
  // which prompted us to make H* = [X,F] rather than H* = [F, X], like in Eidsvik.

  sp_mat SigmaFEandEtaInv = createBlockMatrix(FEinvAndKinvMatrixList) ;
  m_SigmaFEandEtaInvBlockIndices = extractBlockIndices(SigmaFEandEtaInv) ;

  uint basicIndex = 0 ;
  for (uint i = 0 ; i < m_SigmaFEandEtaInvBlockIndices.size() - 1 ; i++) {
    uint blockIndex = m_SigmaFEandEtaInvBlockIndices.at(i) ;
    uint nextBlockIndex = m_SigmaFEandEtaInvBlockIndices.at(i+1) ;
    uint numRows = nextBlockIndex - blockIndex ;
    umat blockPos = join_rows(rep(regspace<uvec>(0, numRows - 1), numRows),
                              rep_each(regspace<uvec>(0, numRows - 1), numRows)) + basicIndex ;

    m_SigmaPos = join_cols(m_SigmaPos, blockPos) ;
    basicIndex += numRows ;
  }
  // uvec keepIndices = find(m_SigmaPos.col(0) >= m_SigmaPos.col(1)) ;
  // m_SigmaPos = m_SigmaPos.rows(keepIndices) ;
  return SigmaFEandEtaInv ;
}

arma::sp_mat AugTree::UpdateSigmaBetaEtaInvMat(Rcpp::Function buildSparse) {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  mat FEinvMatrix = pow(m_fixedEffSD, -2) * eye<mat>(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;
  uint accSize = 0 ;

  vec concatenatedValues(m_SigmaPos.n_rows) ;
  concatenatedValues.subvec(0, size(FEinvMatrix.diag())) = FEinvMatrix.diag() ;
  uint index = FEinvMatrix.n_rows ;
  uint newIndex ;
  for (uint i = 1; i < FEinvAndKinvMatrixList.size(); i++) {
    newIndex = index + FEinvAndKinvMatrixList.at(i)->size() ;
    concatenatedValues.subvec(index,  newIndex - 1) = vectorise(*FEinvAndKinvMatrixList.at(i)) ;
    index = newIndex ;
  }
  uvec dims(2) ;
  dims.fill(GetNumKnots() + m_dataset.covariateValues.n_cols + 1) ;
  sp_mat FEinvAndKinvMatrices = Rcpp::as<sp_mat>(buildSparse(m_SigmaPos, concatenatedValues, dims, false)) ;

  return FEinvAndKinvMatrices ;
}

void AugTree::createHmatrixPos() {

  int numObs = m_dataset.spatialCoords.n_rows ;
  std::vector<std::vector<double *>> memAddressesVec(m_numKnots) ;

  std::vector<uvec> FmatNodeOrder(m_M) ;
  uvec FmatObsOrder(numObs, fill::zeros) ;
  uint rowIndex = 0 ;
  uint colIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

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

    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.subvec(rowIndex, size(observationIndices)) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      uvec brickAncestors(brickList.size()) ;
      for (uint i = 0 ; i < brickAncestors.size() ; i++) {
        brickAncestors.at(i) = brickList.at(i)->GetNodeId() ;
      }
      uint index = 0 ;
      for (uint nodeIndex = 0; nodeIndex < m_vertexVector.size() ; nodeIndex++) {
        int nodePos = GetNodePos(nodeIndex) ;
        if (index < brickAncestors.size()) {
          if (nodeIndex == brickAncestors.at(index)) {
            uvec rowIndices = rep(regspace<uvec>(0, nodeToProcess->GetB(index).n_rows - 1),
                                  nodeToProcess->GetB(index).n_cols) + rowIndex ;
            uvec colIndices = rep_each(regspace<uvec>(0, nodeToProcess->GetB(index).n_cols - 1),
                                       nodeToProcess->GetB(index).n_rows) + colIndex ;

            umat positions = join_rows(rowIndices, colIndices) ;
            m_HmatPos = join_cols(m_HmatPos, positions) ; // Slow, but this is not done many times.

            index += 1 ;
          }
        }
        colIndex += m_vertexVector.at(nodePos)->GetNumKnots() ;
      }
      colIndex = 0 ;
      rowIndex += observationIndices.size() ; // The B matrices should have as many rows as observations in the node...
    }
  }

  m_obsOrderForFmat = FmatObsOrder ;

  mat transIncrementedCovar = trans(join_rows(ones<vec>(numObs), m_dataset.covariateValues)) ;
  m_incrementedCovarReshuffled = trans(transIncrementedCovar.cols(m_obsOrderForFmat)) ;
  m_HmatPos.col(1) += m_dataset.covariateValues.n_cols + 1 ;

  umat covPos = join_rows(rep(regspace<uvec>(0, m_incrementedCovarReshuffled.n_rows - 1), m_incrementedCovarReshuffled.n_cols),
                          rep_each(regspace<uvec>(0, m_incrementedCovarReshuffled.n_cols - 1), m_incrementedCovarReshuffled.n_rows)) ;
  m_HmatPos = join_cols(covPos, m_HmatPos) ;
}

void AugTree::updateHmatrix(Rcpp::Function sparseMatrixConstructFun) {
  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

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

  vec concatenatedValues = vectorise(m_incrementedCovarReshuffled) ;

  for (auto & tipNode : tipNodes) {
    for (auto & Bmat : tipNode->GetWlist()) {
      concatenatedValues = join_cols(concatenatedValues, vectorise(Bmat)) ;
    }
  }

  uvec dims(2) ;
  dims.at(0) = m_dataset.responseValues.size() ;
  dims.at(1) = GetNumKnots() + m_dataset.covariateValues.n_cols + 1 ;

  m_Hmat = Rcpp::as<sp_mat>(sparseMatrixConstructFun(m_HmatPos, concatenatedValues, dims, false)) ;
}

arma::sp_mat AugTree::updateHmatrixPred(Rcpp::Function sparseMatrixConstructFun) {
  vec concatenatedValues = vectorise(m_incrementedCovarPredictReshuffled) ;
  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetPredIndices().size() > 0) {
      for (auto & Umat : tipNode->GetUmatList()) {
        concatenatedValues = join_cols(concatenatedValues, vectorise(Umat)) ;
      }
    }
  }
  uvec dims(2) ;
  dims.at(0) = m_predictData.responseValues.size() ;
  dims.at(1) = m_Hmat.n_cols ;

  return Rcpp::as<sp_mat>(sparseMatrixConstructFun(m_HmatPredPos, concatenatedValues, dims, false)) ;
}

void AugTree::createHmatrixPredPos() {

  uint numObs = m_predictData.spatialCoords.n_rows ;

  std::vector<uvec> FmatNodeOrder(m_M) ;
  uvec FmatObsOrder(numObs, fill::zeros) ;
  uint rowIndex = 0 ;
  uint colIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

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

  m_obsOrderForHpredMat = uvec(m_predictData.timeCoords.size(), fill::zeros) ;
  uint elementIndex = 0 ;
  for (auto & i : tipNodes) {
    uint numPredObsInNode = i->GetPredIndices().size() ;
    if (numPredObsInNode > 0) {
      m_obsOrderForHpredMat.subvec(elementIndex, size(i->GetPredIndices())) = i->GetPredIndices() ;
      elementIndex += i->GetPredIndices().size() ;
    }
  }

  for (auto & nodeToProcess : tipNodes) {
    uvec observationIndices = nodeToProcess->GetPredIndices() ;
    std::vector<TreeNode *> ancestorsList = nodeToProcess->getAncestors() ;

    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.subvec(rowIndex, size(observationIndices)) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      uvec brickAncestors(brickList.size()) ;
      for (uint i = 0 ; i < brickAncestors.size() ; i++) {
        brickAncestors.at(i) = brickList.at(i)->GetNodeId() ;
      }
      uint index = 0 ;

      for (uint nodeIndex = 0; nodeIndex < m_vertexVector.size() ; nodeIndex++) {
        int nodePos = GetNodePos(nodeIndex) ;
        if (index < brickAncestors.size()) {
          if (nodeIndex == brickAncestors.at(index)) {
            uvec rowIndices = rep(regspace<uvec>(0, observationIndices.size() - 1),
                                  ancestorsList.at(index)->GetNumKnots()) + rowIndex ;
            uvec colIndices = rep_each(regspace<uvec>(0, ancestorsList.at(index)->GetNumKnots() - 1),
                                       observationIndices.size()) + colIndex ;
            umat positions = join_rows(rowIndices, colIndices) ;
            m_HmatPredPos = join_cols(m_HmatPredPos, positions) ; // Slow, but this is not done many times.

            index += 1 ;
          }
        }
        colIndex += m_vertexVector.at(nodePos)->GetNumKnots() ;
      }
      colIndex = 0 ;
      rowIndex += observationIndices.size() ;
    }
  }
  mat transIncrementedCovar = trans(join_rows(ones<vec>(numObs), m_predictData.covariateValues)) ;
  m_incrementedCovarPredictReshuffled = trans(transIncrementedCovar.cols(FmatObsOrder)) ;

  m_HmatPredPos.col(1) += m_predictData.covariateValues.n_cols + 1 ;
  umat covPos = join_rows(rep(regspace<uvec>(0, m_incrementedCovarPredictReshuffled.n_rows - 1),
                              m_incrementedCovarPredictReshuffled.n_cols),
                          rep_each(regspace<uvec>(0, m_incrementedCovarPredictReshuffled.n_cols - 1),
                                   m_incrementedCovarPredictReshuffled.n_rows)) ;
  m_HmatPredPos = join_cols(covPos, m_HmatPredPos) ;
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

  std::vector<std::pair<double, GammaHyperParas>> priorCombinations ;

  priorCombinations.push_back(std::make_pair(m_MRAcovParasSpace.m_rho, m_maternParasGammaAlphaBetaSpace.m_rho)) ;
  priorCombinations.push_back(std::make_pair(m_MRAcovParasSpace.m_smoothness, m_maternParasGammaAlphaBetaSpace.m_smoothness)) ;

  priorCombinations.push_back(std::make_pair(m_MRAcovParasTime.m_rho, m_maternParasGammaAlphaBetaTime.m_rho)) ;
  priorCombinations.push_back(std::make_pair(m_MRAcovParasTime.m_smoothness, m_maternParasGammaAlphaBetaTime.m_smoothness)) ;

  priorCombinations.push_back(std::make_pair(m_spacetimeScaling, m_maternSpacetimeScalingGammaAlphaBeta)) ;

  priorCombinations.push_back(std::make_pair(m_fixedEffSD, m_fixedEffGammaAlphaBeta)) ;
  priorCombinations.push_back(std::make_pair(m_errorSD, m_errorGammaAlphaBeta)) ;

  double logPrior = 0 ;

  for (auto & i : priorCombinations) {
    logPrior += (i.second.m_alpha - 1) * log(i.first) - i.second.m_beta * i.first ;
  }

  m_logPrior = logPrior ;
}

void AugTree::ComputeLogFCandLogCDandDataLL(Rcpp::Function gradCholeskiFun, Rcpp::Function sparseMatrixConstructFun,
                                            Rcpp::Function sparseDeterminantFun) {

  int n = m_dataset.responseValues.size() ;
  // cout << "Creating matrix of Ks... \n" ;
  sp_mat SigmaFEandEtaInv ;
  if (m_SigmaPos.n_rows == 0) {
    SigmaFEandEtaInv = CreateSigmaBetaEtaInvMat() ;
  } else {
    SigmaFEandEtaInv = UpdateSigmaBetaEtaInvMat(sparseMatrixConstructFun) ;
  }
  double logDetSigmaKFEinv = logDetBlockMatrix(SigmaFEandEtaInv, m_SigmaFEandEtaInvBlockIndices) ;
  // cout << "Done... \n" ;

  if (m_recomputeMRAlogLik) {
    // cout << "Obtaining H matrix... \n" ;
    if (m_HmatPos.size() == 0) {
      createHmatrixPos() ;
    }
    updateHmatrix(sparseMatrixConstructFun) ;

    // cout << "Done... \n" ;
  }

  sp_mat secondTerm = std::pow(m_errorSD, -2) * trans(m_Hmat) * m_Hmat ;
  // cout << "Obtaining Q... \n" ;
  m_FullCondPrecision = SigmaFEandEtaInv + secondTerm ;

  // cout << "Done... \n" ;

  vec responsesReshuffled = m_dataset.responseValues.elem(m_obsOrderForFmat) ;

  // cout << "Computing FC mean... \n" ;

  // sp_mat hessianMat = SigmaFEandEtaInv + secondTerm ;
  mat scaledResponse = std::pow(m_errorSD, -2) * trans(responsesReshuffled) * m_Hmat ;

  Rcpp::NumericVector updatedMean = gradCholeskiFun(m_FullCondPrecision, scaledResponse) ;

  m_Vstar = updatedMean ; // Assuming there will be an implicit conversion to vec type.
  m_FullCondMean = m_Vstar ;

  vec fixedEffMeans = m_Vstar.head(m_fixedEffParameters.size()) ;
  SetFixedEffParameters(fixedEffMeans) ;

  // double logDetQmat = logDeterminantQmat(sparseMatrixConstructFun) ;
  double logDetQmat = Rcpp::as<double>(sparseDeterminantFun(m_FullCondPrecision)) ;

  m_logFullCond = 0.5 * logDetQmat ; // Since we arbitrarily evaluate always at the full-conditional mean, the exponential part of the distribution reduces to 0.

  // Computing p(v* | Psi)
  mat logLikTerm = -0.5 * trans(m_Vstar) * SigmaFEandEtaInv * m_Vstar ;

  m_logCondDist = logLikTerm(0,0) + 0.5 * logDetSigmaKFEinv ;

  // Computing p(y | v*, Psi)
  double errorLogDet = -2 * n * log(m_errorSD) ;
  vec recenteredY = responsesReshuffled - m_Hmat * m_Vstar ;
  vec globalLogLikExp = -0.5 * std::pow(m_errorSD, -2) * trans(recenteredY) * recenteredY ;
  m_globalLogLik = 0.5 * errorLogDet + globalLogLikExp(0) ;
}

void AugTree::ComputeLogJointPsiMarginal(Rcpp::Function gradCholeskiFun, Rcpp::Function sparseMatConstructFun,
                                         Rcpp::Function sparseDeterminantFun) {

  ComputeLogPriors() ;
  // cout << "Computing Wmats... \n" ;
  if (m_recomputeMRAlogLik) {
    m_MRAcovParasSpace.print("Space parameters:") ;
    m_MRAcovParasTime.print("Time parameters:") ;
    printf("Scale parameter: %.4e \n", m_spacetimeScaling) ;
    fflush(stdout); // Will now print everything in the stdout buffer
    computeWmats() ; // This will produce the K matrices required. NOTE: ADD CHECK THAT ENSURES THAT THE MRA LIK. IS ONLY RE-COMPUTED WHEN THE MRA COV. PARAMETERS CHANGE.
    TreeNode * arbitraryTipNode = GetTipNodes().at(0) ;
    std::vector<TreeNode *> ancestorsList = arbitraryTipNode->getAncestors() ;
    fflush(stdout) ;
    cout << "\n W_j1^0 \n\n" ;
    printf("Dimensions: %i %i \n\n", ancestorsList.at(1)->GetWlist().at(0).n_rows, ancestorsList.at(1)->GetWlist().at(1).n_cols) ;
    std::cout << ancestorsList.at(1)->GetWlist().at(0)(0,0,size(3,3)) << "\n\n" ;
    printf("Ancestors node ID: %i \n", ancestorsList.at(1)->GetNodeId()) ;
    // std::cout << "These are the knots at resolution 0: \n" ;
    // std::cout << ancestorsList.at(0)->GetKnotsCoor().spatialCoords(0,0,size(10,2)) << "\n\n";
    // std::cout << ancestorsList.at(1)->GetKnotsCoor().timeCoords.subvec(0,9) << "\n\n";
    // std::cout << "These are the knots at resolution 1: \n" ;
    // std::cout << ancestorsList.at(0)->GetKnotsCoor().spatialCoords(0,0,size(20,2)) << "\n\n";
    // std::cout << ancestorsList.at(1)->GetKnotsCoor().timeCoords.subvec(0,19) << "\n\n";
  }

  ComputeLogFCandLogCDandDataLL(gradCholeskiFun, sparseMatConstructFun, sparseDeterminantFun) ;

  printf("Observations log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
  m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  m_logJointPsiMarginal = m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond ;
  printf("Joint value: %.4e \n \n", m_logJointPsiMarginal) ;
}

// This inversion is based on recursive partitioning of the Q matrix. It is based on the observation that it is
// possible to form block-diagonal matrices on the diagonal which can be easily inverted.
// The challenge in inverting Qmat was the very heavy memory burden.
// This function involves much smaller matrices, which will make the operations easier to handle.
// double AugTree::logDeterminantQmat(Rcpp::Function funToConstructSparse) {
//
//   if (m_DmatrixBlockIndices.size() == 0) {
//     m_DmatrixBlockIndices = extractBlockIndicesFromLowerRight(m_FullCondPrecision) ;
//
//     uint basicIndex = 0 ;
//     for (uint i = 0 ; i < m_DmatrixBlockIndices.size() - 1 ; i++) {
//       uint blockIndex = m_DmatrixBlockIndices.at(i) ;
//       uint nextBlockIndex = m_DmatrixBlockIndices.at(i+1) ;
//       uint numRows = nextBlockIndex - blockIndex ;
//       umat blockPos = join_rows(rep(regspace<uvec>(0, numRows - 1), numRows),
//                                 rep_each(regspace<uvec>(0, numRows - 1), numRows)) + basicIndex ;
//       m_DinFCmatPos = join_cols(m_DinFCmatPos, blockPos) ;
//       basicIndex += numRows ;
//     }
//   }
//
//   int numRowsD = m_FullCondPrecision.n_rows - m_DmatrixBlockIndices(0) ;
//   int numRowsA = m_DmatrixBlockIndices(0) ;
//
//   int shift = m_DmatrixBlockIndices.at(0) ;
//
//   uvec shiftedBlockIndices = m_DmatrixBlockIndices - shift ;
//
//   vec concatenatedValues ;
//
//   for (uint i = 0 ; i < m_DmatrixBlockIndices.size() - 1 ; i++) {
//     uint index = m_DmatrixBlockIndices.at(i) ;
//     uint numRows = m_DmatrixBlockIndices.at(i + 1) - index ;
//     vec vectorisedInverse = vectorise(inv_sympd(mat(m_FullCondPrecision(index, index, size(numRows, numRows))))) ;
//     concatenatedValues = join_cols(concatenatedValues, vectorisedInverse) ;
//   }
//
//   sp_mat Dinv = Rcpp::as<sp_mat>(funToConstructSparse(m_DinFCmatPos, concatenatedValues)) ;
//
//   double logDeterminantD = logDetBlockMatrix(m_FullCondPrecision(m_DmatrixBlockIndices.at(0), m_DmatrixBlockIndices.at(0), size(numRowsD, numRowsD)), shiftedBlockIndices) ;
//
//   sp_mat Amatrix = m_FullCondPrecision(0, 0, size(numRowsA, numRowsA)) ;
//
//   double logDeterminantComposite, sign1 ;
//   // uint AmatrixSize = DmatrixBlockIndices.at(0) ;
//
//   sp_mat Bmatrix = m_FullCondPrecision(0, m_DmatrixBlockIndices.at(0), size(numRowsA, numRowsD)) ;
//   // The next few lines ensure that the matrix whose determinant needs to be computed is
//   // really symmetric. Else computational zeros will ruin the symmetry.
//   sp_mat compositeMat = Amatrix - Bmatrix * Dinv * trans(Bmatrix) ;
//
//   log_det(logDeterminantComposite, sign1, mat(compositeMat)) ;
//
//   if (sign1 < 0) {
//     throw Rcpp::exception("Error in logDeterminantQmat! sign1 should be positive. \n") ;
//   } // The determinant for the composite must be positive, because the determinant for D is positive. If it were negative, the determinant for the Q matrix would be negative, which is not allowed since it's a covariance matrix.
//
//   double logDeterminant = logDeterminantD + logDeterminantComposite ;
//   // cout << "Leaving logDeterminantQmat... \n" ;
//   return logDeterminant ;
// }

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

sp_mat AugTree::ComputeHpred(const mat & spCoords, const vec & time, const mat & covariateMatrix, Rcpp::Function sparseMatrixConstructFun) {

  // SetPredictData(spCoords, time, covariateMatrix) ;
  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

  for (auto & i : tipNodes) {
    // i->SetPredictLocations(m_predictData) ;
    uvec predictionsInLeaf = i->GetPredIndices() ;
    if (predictionsInLeaf.size() > 0) {
      i->computeUpred(m_MRAcovParasSpace, m_MRAcovParasTime, m_spacetimeScaling, m_predictData, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;
    }
  }
  if (m_HmatPredPos.n_rows == 0) {
    createHmatrixPredPos() ;
  }

  return updateHmatrixPred(sparseMatrixConstructFun) ;
}

arma::vec AugTree::ComputeEvar(const arma::sp_mat & HmatPred, Rcpp::Function sparseSolveFun, const int batchSize) {

  double errorVar = std::pow(m_errorSD, 2) ;
  vec EvarValues(HmatPred.n_rows, fill::zeros) ;
  int obsIndex = 0 ;

  while (obsIndex < HmatPred.n_rows) {

    int newObsIndex = std::min(obsIndex + batchSize - 1, int(HmatPred.n_rows) - 1) ;
    sp_mat bVec = HmatPred.rows(obsIndex, newObsIndex) ;
    sp_mat bVecTrans = trans(HmatPred.rows(obsIndex, newObsIndex)) ;
    sp_mat meanValue = bVecTrans % Rcpp::as<sp_mat>(sparseSolveFun(m_FullCondPrecision, bVecTrans)) ;
    EvarValues.subvec(obsIndex, newObsIndex) = trans(sum(meanValue,0)) + errorVar ;

    obsIndex += batchSize ;
  }
  return EvarValues ;
}

void AugTree::SetMRAcovParas(const Rcpp::List & MRAcovParas) {
  List SpaceParas = Rcpp::as<List>(MRAcovParas["space"]) ;
  List TimeParas = Rcpp::as<List>(MRAcovParas["time"]) ;
  double scalePara = Rcpp::as<double>(MRAcovParas["scale"]) ;

  double rhoSpace = Rcpp::as<double>(SpaceParas["rho"]) ;
  double smoothnessSpace = Rcpp::as<double>(SpaceParas["smoothness"]) ;

  double rhoTime = Rcpp::as<double>(TimeParas["rho"]) ;
  double smoothnessTime = Rcpp::as<double>(TimeParas["smoothness"]) ;

  maternVec MRAcovParasSpace(rhoSpace, smoothnessSpace, 1) ;
  maternVec MRAcovParasTime(rhoTime, smoothnessTime, 1) ;

  bool test = (fabs(m_spacetimeScaling - scalePara) < epsilon) && (m_MRAcovParasSpace == MRAcovParasSpace) && (m_MRAcovParasTime == MRAcovParasTime) ;
  // I overloaded == to have it use the epsilon threshold as well.

  if (test) {
    m_recomputeMRAlogLik = false ;
  } else {
    m_recomputeMRAlogLik = true ;
  }
  m_MRAcovParasSpace = MRAcovParasSpace ;
  m_MRAcovParasTime = MRAcovParasTime ;
  m_spacetimeScaling = scalePara ;
}

void AugTree::SetMRAcovParasGammaAlphaBeta(const Rcpp::List & MRAcovParasList) {
  Rcpp::List spaceParas = Rcpp::as<List>(MRAcovParasList["space"]) ;
  Rcpp::List timeParas = Rcpp::as<List>(MRAcovParasList["time"]) ;
  vec scaling = Rcpp::as<vec>(MRAcovParasList["scale"]) ;
  m_maternSpacetimeScalingGammaAlphaBeta = scaling ;
  m_maternParasGammaAlphaBetaSpace = maternGammaPriorParasWithoutScale(GammaHyperParas(Rcpp::as<vec>(spaceParas["rho"])),
                                            GammaHyperParas(Rcpp::as<vec>(spaceParas["smoothness"]))) ;
  m_maternParasGammaAlphaBetaTime = maternGammaPriorParasWithoutScale(GammaHyperParas(Rcpp::as<vec>(timeParas["rho"])),
                                            GammaHyperParas(Rcpp::as<vec>(timeParas["smoothness"]))) ;
}
