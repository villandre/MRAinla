#include <math.h>
// #include <omp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"

using namespace Rcpp ;
using namespace MRAinla ;
using namespace Eigen ;

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

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, vec & timeRange, vec & observations, mat & obsSp, vec & obsTime, mat & predCovariates, mat & predSp, vec & predTime, uint & minObsForTimeSplit, unsigned long int & seed, mat & covariates, const bool splitTime, const unsigned int numKnotsRes0, const unsigned int J)
  : m_M(M)
{
  m_dataset = inputdata(observations, obsSp, obsTime, covariates) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  m_randomNumGenerator = gsl_rng_alloc(gsl_rng_taus) ;
  SetPredictData(predSp, predTime, predCovariates) ;

  gsl_rng_set(m_randomNumGenerator, seed) ;
  // m_fixedEffParameters.resize(m_dataset.covariateValues.cols() + 1) ; // An extra 1 is added to take into account the intercept.
  m_fixedEffParameters = Eigen::VectorXd::Zero(m_dataset.covariateValues.cols() + 1);

  BuildTree(minObsForTimeSplit, splitTime, numKnotsRes0, J) ;

  m_numKnots = 0 ;

  for (uint i = 0 ; i < m_vertexVector.size() ; i++) {
    m_numKnots += m_vertexVector.at(i)->GetKnotsCoor().timeCoords.size() ;
  }
  m_FullCondMean = Eigen::VectorXd::Zero(m_numKnots + m_fixedEffParameters.size()) ;
  m_FullCondSDs = Eigen::VectorXd::Zero(m_numKnots + m_fixedEffParameters.size()) ;
  m_Vstar = Eigen::VectorXd::Zero(m_numKnots + m_dataset.covariateValues.cols() + 1) ;
  m_Vstar.segment(0, m_FEmu.size()) = m_FEmu ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

  for (auto & i : tipNodes) {
    i->SetPredictLocations(m_predictData) ;
  }

  createHmatrixPredPos() ;
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
  Eigen::PermutationMatrix<Dynamic, Dynamic> perm(obsForMedian) ;
  uvec childMembership = uvec::Zero(obsForMedian.size()) ;
  std::vector<dimensions> childDimensions ;
  childDimensions.push_back(parent->GetDimensions()) ;

  if (obsForMedian.size() <= 1) {
    throw Rcpp::exception("Cannot split empty region or region with only one observation.\n") ;
  }

  uint newChildIndex = 1 ;

  // Handling longitude
  uvec currentChildIndices = uvec::LinSpaced(newChildIndex, 0, newChildIndex - 1) ;
  // for (auto & j : currentChildIndices) {
  for (uint j = 0 ; j < currentChildIndices.size(); j++) {
    uint currentChildIndex = currentChildIndices(j) ;
    Eigen::Matrix<bool, Dynamic, 1> compareVec = childMembership == currentChildIndex ;
    uvec elementsInChild = find(compareVec) ;
    vec column = m_dataset.spatialCoords.col(0) ;
    vec subColumn = perm * column ;
    vec subVector(elementsInChild.size()) ;
    for (uint innerIndex = 0; innerIndex < elementsInChild.size(); innerIndex++) {
      subVector(innerIndex) = subColumn(elementsInChild(innerIndex)) ;
    }
    double colMedian = median(subVector) ;
    vec updatedLongitude{childDimensions.at(currentChildIndices(j)).longitude(0), colMedian} ;
    vec newChildLongitude{colMedian, childDimensions.at(currentChildIndices(j)).longitude(1)} ;
    dimensions newDimensions = childDimensions.at(currentChildIndices(j)) ;
    newDimensions.longitude = newChildLongitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(j).longitude = updatedLongitude ;
    uvec greaterElements = find(elem(subColumn, elementsInChild) > colMedian) ; // In deriveObsInNode, the checks are <=. It follows that observations on the right and upper boundaries of a zone are included.
    PermutationMatrix<Dynamic, Dynamic> perm(greaterElements) ;
    uvec updateIndices = perm * elementsInChild ;
    for (uint innerIndex = 0; innerIndex < updateIndices.size(); innerIndex++) {
      childMembership(updateIndices(innerIndex)) = newChildIndex ;
    }
    newChildIndex += 1 ;
  }

  // Handling latitude
  currentChildIndices = uvec::LinSpaced(newChildIndex, 0, newChildIndex - 1) ;
  for (uint j = 0 ; j < currentChildIndices.size(); j++) {
    vec column = m_dataset.spatialCoords.col(1) ;
    vec subColumn = perm * column ;
    uvec elementsInChild = find(childMembership == currentChildIndices(j)) ;
    double colMedian = median(elem(subColumn, elementsInChild)) ;
    vec updatedLatitude{childDimensions.at(currentChildIndices(j)).latitude(0), colMedian} ;
    vec newChildLatitude{colMedian, childDimensions.at(currentChildIndices(j)).latitude(1)} ;
    dimensions newDimensions = childDimensions.at(currentChildIndices(j)) ;
    newDimensions.latitude = newChildLatitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(currentChildIndices(j)).latitude = updatedLatitude ;
    uvec greaterElements = find(elem(subColumn, elementsInChild) > colMedian) ;
    uvec updateIndices = elem(elementsInChild, greaterElements) ;
    for (uint innerIndex = 0; innerIndex < updateIndices.size(); innerIndex++) {
      childMembership(updateIndices(innerIndex)) = newChildIndex ;
    }
    newChildIndex += 1 ;
  }
  if (splitTime) {
    // Handling time
    currentChildIndices = uvec::LinSpaced(newChildIndex, 0, newChildIndex - 1) ;
    for (uint j = 0 ; j < currentChildIndices.size(); j++) {
      vec time = perm * m_dataset.timeCoords ;
      uvec elementsInChild = find(childMembership == currentChildIndices(j)) ;
      vec elementsForMedian = elem(time, elementsInChild) ;
      uvec uniqueTimeValues = find_unique(elementsForMedian) ;

      if (uniqueTimeValues.size() > 1) {
        double colMedian = median(elementsForMedian) ;
        vec updatedTime{childDimensions.at(currentChildIndices(j)).time(0), colMedian} ;
        vec newChildTime{colMedian, childDimensions.at(currentChildIndices(j)).time(1)} ;
        dimensions newDimensions = childDimensions.at(currentChildIndices(j)) ;
        newDimensions.time = newChildTime ;
        childDimensions.push_back(newDimensions) ;
        childDimensions.at(currentChildIndices(j)).time = updatedTime ;
        uvec greaterElements = find(elem(time, elementsInChild) > colMedian) ;
        Eigen::PermutationMatrix<Dynamic, Dynamic> permGreaterElements(greaterElements) ;
        uvec updateIndices = permGreaterElements * greaterElements ;
        for (uint innerIndex = 0 ; innerIndex < updateIndices.size() ; innerIndex++) {
          childMembership(updateIndices(innerIndex)) = newChildIndex ;
        }
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

void AugTree::ComputeMRAlogLikAlt(const bool WmatsAvailable)
{
  if (!WmatsAvailable) {
    computeWmats() ;
  }
  double MRAlogLik = 0 ;
  int currentIndex = m_fixedEffParameters.size() ;

  for (auto & i : m_vertexVector) {
    int lastIndex = currentIndex + i->GetNumKnots() - 1 ;
    vec etas = m_Vstar.segment(currentIndex, lastIndex - currentIndex + 1) ;

    double logDeterminantValue = 0 ;
    double sign = 0 ;
    logDeterminantValue = i->GetKmatrix().colPivHouseholderQr().logAbsDeterminant() ;

    mat exponentTerm = -0.5 * etas.transpose() * i->GetKmatrixInverse() * etas ;
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

    // #pragma omp parallel for
    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++) {
      (*it)->DeriveAtilde() ;
    }
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

sp_mat AugTree::CreateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  mat FEinvMatrix = pow(m_fixedEffSD, -2) * mat::Identity(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;

  // This is because we want the Q matrix to have a block-diagonal component in the lower-right corner,
  // which prompted us to make H* = [X,F] rather than H* = [F, X], like in Eidsvik.

  sp_mat SigmaFEandEtaInv = createBlockMatrix(FEinvAndKinvMatrixList) ;
  m_SigmaFEandEtaInvBlockIndices = extractBlockIndices(SigmaFEandEtaInv) ;

  uint basicIndex = 0 ;
  for (uint i = 0 ; i < m_SigmaFEandEtaInvBlockIndices.size() - 1 ; i++) {
    uint blockIndex = m_SigmaFEandEtaInvBlockIndices(i) ;
    uint nextBlockIndex = m_SigmaFEandEtaInvBlockIndices(i+1) ;
    uint numRows = nextBlockIndex - blockIndex ;
    umat blockPos = join_rows(rep(uvec::LinSpaced(numRows, 0, numRows - 1), numRows),
                              rep_each(uvec::LinSpaced(numRows, 0, numRows - 1), numRows)) + basicIndex ;

    m_SigmaPos = join_cols(m_SigmaPos, blockPos) ;
    basicIndex += numRows ;
  }

  return SigmaFEandEtaInv ;
}

sp_mat AugTree::UpdateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;

  mat FEinvMatrix = pow(m_fixedEffSD, -2) * mat::Identity(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;

  uint numElements = 0 ;
  for (auto & i : FEinvAndKinvMatrixList) { // This would not need to be done everytime an update is processed, as numElements doesn't change, but it's probably very fast, so it shouldn't matter.
    numElements += i->size() ;
  }

  vec concatenatedValues(numElements) ;
  concatenatedValues.segment(0, m_fixedEffParameters.size()) = FEinvMatrix.diagonal() ;
  uint index = FEinvMatrix.rows() ;

  for (uint i = 1; i < FEinvAndKinvMatrixList.size(); i++) {
    concatenatedValues.segment(index,  FEinvAndKinvMatrixList.at(i)->size()) = Map<vec>((*FEinvAndKinvMatrixList.at(i)).data(), (*FEinvAndKinvMatrixList.at(i)).cols() * (*FEinvAndKinvMatrixList.at(i)).rows()) ;
    index += FEinvAndKinvMatrixList.at(i)->size() ;
  }

  int loopIndex = 0 ;
  sp_mat FEinvAndKinvMatrices(numElements, numElements);
  for (int k = 0; k < FEinvAndKinvMatrices.outerSize(); ++k) {
    for (SparseMatrix<double>::InnerIterator it(FEinvAndKinvMatrices, k); it; ++it)
    {
      it.valueRef() = concatenatedValues(loopIndex) ;
      loopIndex += 1 ;
    }
  }

  return FEinvAndKinvMatrices ;
}

void AugTree::createHmatrixPos() {

  int numObs = m_dataset.spatialCoords.rows() ;
  std::vector<std::vector<double *>> memAddressesVec(m_numKnots) ;

  std::vector<uvec> FmatNodeOrder(m_M) ;
  uvec FmatObsOrder = uvec::Zero(numObs) ;
  int rowIndex = 0 ;
  int colIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
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
      bool test = (firstAncestorIds(i) == secondAncestorIds(i)) ;
      if (!test) {
        firstSmallerThanSecond = (firstAncestorIds(i) < secondAncestorIds(i)) ;
        break ;
      }
    }
    return firstSmallerThanSecond ;
  }) ; // This is supposed to reorder the tip nodes in such a way that the F matrix has contiguous blocks.

  for (auto & nodeToProcess : tipNodes) {
    uvec observationIndices = nodeToProcess->GetObsInNode() ;

    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.segment(rowIndex, observationIndices.size()) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      uvec brickAncestors(brickList.size()) ;
      for (uint i = 0 ; i < brickAncestors.size() ; i++) {
        brickAncestors(i) = brickList.at(i)->GetNodeId() ;
      }
      uint index = 0 ;
      for (uint nodeIndex = 0; nodeIndex < m_vertexVector.size() ; nodeIndex++) {
        int nodePos = GetNodePos(nodeIndex) ;
        if (index < brickAncestors.size()) {
          if (nodeIndex == brickAncestors(index)) {
            uvec rowIndices = rep(uvec::LinSpaced(nodeToProcess->GetB(index).rows(), 0, nodeToProcess->GetB(index).rows() - 1),
                                  nodeToProcess->GetB(index).cols()) + rowIndex ;
            uvec colIndices = rep_each(uvec::LinSpaced(nodeToProcess->GetB(index).cols(), 0, nodeToProcess->GetB(index).cols() - 1),
                                       nodeToProcess->GetB(index).rows()) + colIndex ;

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

  Eigen::PermutationMatrix<Dynamic, Dynamic> perm(m_obsOrderForFmat) ;
  mat transIncrementedCovar = (join_rows(vec::Ones(numObs), m_dataset.covariateValues)).tranpose() ;
  m_incrementedCovarReshuffled = (transIncrementedCovar * perm).transpose() ;
  m_HmatPos.col(1) += m_dataset.covariateValues.cols() + 1 ;

  umat covPos = join_rows(rep(uvec::LinSpaced(m_incrementedCovarReshuffled.rows(), 0, m_incrementedCovarReshuffled.rows() - 1), m_incrementedCovarReshuffled.cols()),
                          rep_each(uvec::LinSpaced(m_incrementedCovarReshuffled.cols(), 0, m_incrementedCovarReshuffled.cols() - 1), m_incrementedCovarReshuffled.rows())) ;
  m_HmatPos = join_cols(covPos, m_HmatPos) ;
}

void AugTree::updateHmatrix() {
  int numElementsInH = m_incrementedCovarReshuffled.size() ;
  for (auto & tipNode : GetTipNodes()) {
    for (auto & Bmat : tipNode->GetWlist()) {
      numElementsInH += Bmat.size() ;
    }
  }
  vec concatenatedValues(numElementsInH) ;
  concatenatedValues.segment(0, m_incrementedCovarReshuffled.size()) = Map<vec>(m_incrementedCovarReshuffled, m_incrementedCovarReshuffled.size()) ;
  int index = m_incrementedCovarReshuffled.size() ;
  for (auto & tipNode : GetTipNodes()) {
    for (auto & Bmat : tipNode->GetWlist()) {
      concatenatedValues.segment(index, Bmat.size()) = Map<vec>(Bmat, Bmat.size()) ;
      index += Bmat.size() ;
    }
  }

  int loopIndex = 0 ;
  for (int k=0; k<m_Hmat.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_Hmat, k); it; ++it)
    {
      it.valueRef() = concatenatedValues(loopIndex) ;
      loopIndex += 1 ;
    }
  }
}

void AugTree::updateHmatrixPred() {
  vec concatenatedValues = vectorise(m_incrementedCovarPredictReshuffled) ;
  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetPredIndices().size() > 0) {
      for (auto & Umat : tipNode->GetUmatList()) {
        concatenatedValues = join_cols(concatenatedValues, vectorise(Umat)) ;
      }
    }
  }

  int loopIndex = 0 ;
  for (int k = 0; k < m_HmatPred.outerSize(); ++k) {
    for (SparseMatrix<double>::InnerIterator it(m_HmatPred, k); it; ++it)
    {
      it.valueRef() = concatenatedValues(loopIndex) ;
      loopIndex += 1 ;
    }
  }
}

void AugTree::createHmatrixPredPos() {

  uint numObs = m_predictData.spatialCoords.rows() ;

  std::vector<uvec> FmatNodeOrder(m_M) ;
  uvec FmatObsOrder = uvec::Zero(numObs) ;
  uint rowIndex = 0 ;
  uint colIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
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
      bool test = (firstAncestorIds(i) == secondAncestorIds(i)) ;
      if (!test) {
        firstSmallerThanSecond = (firstAncestorIds(i) < secondAncestorIds(i)) ;
        break ;
      }
    }
    return firstSmallerThanSecond ;
  }) ; // This is supposed to reorder the tip nodes in such a way that the F matrix has contiguous blocks.

  m_obsOrderForHpredMat = uvec::Zero(m_predictData.timeCoords.size()) ;
  uint elementIndex = 0 ;
  for (auto & i : tipNodes) {
    uint numPredObsInNode = i->GetPredIndices().size() ;
    if (numPredObsInNode > 0) {
      m_obsOrderForHpredMat.segment(elementIndex, i->GetPredIndices().size()) = i->GetPredIndices() ;
      elementIndex += i->GetPredIndices().size() ;
    }
  }

  for (auto & nodeToProcess : tipNodes) {
    uvec observationIndices = nodeToProcess->GetPredIndices() ;
    std::vector<TreeNode *> ancestorsList = nodeToProcess->getAncestors() ;

    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.segment(rowIndex, observationIndices.size()) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      uvec brickAncestors(brickList.size()) ;
      for (uint i = 0 ; i < brickAncestors.size() ; i++) {
        brickAncestors(i) = brickList.at(i)->GetNodeId() ;
      }
      uint index = 0 ;

      for (uint nodeIndex = 0; nodeIndex < m_vertexVector.size() ; nodeIndex++) {
        int nodePos = GetNodePos(nodeIndex) ;
        if (index < brickAncestors.size()) {
          if (nodeIndex == brickAncestors(index)) {
            uvec rowIndices = rep(uvec::LinSpaced(observationIndices.size(), 0, observationIndices.size() - 1),
                                  ancestorsList.at(index)->GetNumKnots()) + rowIndex ;
            uvec colIndices = rep_each(uvec::LinSpaced(ancestorsList.at(index)->GetNumKnots(), 0, ancestorsList.at(index)->GetNumKnots() - 1),
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
  Eigen::PermutationMatrix<Dynamic, Dynamic> perm(FmatObsOrder) ;

  mat transIncrementedCovar = (join_rows(vec::Ones(numObs), m_predictData.covariateValues)).transpose() ;
  m_incrementedCovarPredictReshuffled = (transIncrementedCovar * perm).transpose() ;

  m_HmatPredPos.col(1) += m_predictData.covariateValues.cols() + 1 ;
  umat covPos = join_rows(rep(uvec::LinSpaced(m_incrementedCovarPredictReshuffled.rows(), 0, m_incrementedCovarPredictReshuffled.rows() - 1),
                              m_incrementedCovarPredictReshuffled.cols()),
                          rep_each(uvec::LinSpaced(m_incrementedCovarPredictReshuffled.cols(), 0, m_incrementedCovarPredictReshuffled.cols() - 1),
                                   m_incrementedCovarPredictReshuffled.rows())) ;
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

void AugTree::ComputeLogFCandLogCDandDataLL() {

  int n = m_dataset.responseValues.size() ;
  // cout << "Creating matrix of Ks... \n" ;
  sp_mat SigmaFEandEtaInv ;
  if (m_SigmaPos.n_rows == 0) {
    SigmaFEandEtaInv = CreateSigmaBetaEtaInvMat() ;
  } else {
    SigmaFEandEtaInv = UpdateSigmaBetaEtaInvMat() ;
  }
  // double logDetSigmaKFEinv = logDetBlockMatrix(SigmaFEandEtaInv, m_SigmaFEandEtaInvBlockIndices) ;
  Eigen::CholmodSimplicialLDLT<sp_mat> blockChol(SigmaFEandEtaInv) ;
  double logDetSigmaKFEinv = blockChol.logDeterminant() ;
  cout << "Done... \n" ;

  if (m_recomputeMRAlogLik) {
    // cout << "Obtaining H matrix... \n" ;
    if (m_HmatPos.size() == 0) {
      createHmatrixPos() ;
    }
    updateHmatrix() ;

    // cout << "Done... \n" ;
  }

  sp_mat secondTerm = std::pow(m_errorSD, -2) * m_Hmat.transpose() * m_Hmat ;
  vec responsesReshuffled = elem(m_dataset.responseValues, m_obsOrderForFmat) ;

  // sp_mat hessianMat = SigmaFEandEtaInv + secondTerm ;
  mat scaledResponse = std::pow(m_errorSD, -2) * responsesReshuffled.transpose() * m_Hmat ;
  cout << "Obtaining Qchol... \n" ;
  m_FullCondPrecisionChol = Eigen::SimplicialCholesky<sp_mat>::compute(SigmaFEandEtaInv + secondTerm) ;
  cout << "Done! \n" ;
  // m_FullCondSDs = sqrt(m_FullCondPrecisionChol.diag()) ;
  // cout << "Done... \n" ;

  Rcpp::NumericVector updatedMean = m_FullCondPrecisionChol.solve(scaledResponse) ;

  m_Vstar = updatedMean ; // Assuming there will be an implicit conversion to vec type.

  m_FullCondMean = m_Vstar ;

  vec fixedEffMeans = m_Vstar.head(m_fixedEffParameters.size()) ;
  SetFixedEffParameters(fixedEffMeans) ;

  double logDetQmat = m_FullCondPrecisionChol.logDeterminant() ;
  m_logFullCond = 0.5 * logDetQmat ; // Since we arbitrarily evaluate always at the full-conditional mean, the exponential part of the distribution reduces to 0.

  // Computing p(v* | Psi)
  mat logLikTerm = -0.5 * m_Vstar.transpose() * SigmaFEandEtaInv * m_Vstar ;

  m_logCondDist = logLikTerm(0,0) + 0.5 * logDetSigmaKFEinv ;

  // Computing p(y | v*, Psi)
  double errorLogDet = -2 * n * log(m_errorSD) ;
  vec recenteredY = responsesReshuffled - m_Hmat * m_Vstar ;
  vec globalLogLikExp = -0.5 * std::pow(m_errorSD, -2) * recenteredY.transpose() * recenteredY ;
  m_globalLogLik = 0.5 * errorLogDet + globalLogLikExp(0) ;
}

void AugTree::ComputeLogJointPsiMarginal() {

  ComputeLogPriors() ;
  // cout << "Computing Wmats... \n" ;
  if (m_recomputeMRAlogLik) {
    // m_MRAcovParasSpace.print("Space parameters:") ;
    // m_MRAcovParasTime.print("Time parameters:") ;
    // printf("Scale parameter: %.4e \n", m_spacetimeScaling) ;
    // fflush(stdout); // Will now print everything in the stdout buffer
    computeWmats() ; // This will produce the K matrices required. NOTE: ADD CHECK THAT ENSURES THAT THE MRA LIK. IS ONLY RE-COMPUTED WHEN THE MRA COV. PARAMETERS CHANGE.
    // cout << "Done... \n" ;
  }

  ComputeLogFCandLogCDandDataLL() ;

  // printf("Observations log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
  // m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  m_logJointPsiMarginal = m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond ;
  // printf("Joint value: %.4e \n \n", m_logJointPsiMarginal) ;
}

uvec AugTree::extractBlockIndicesFromLowerRight(const sp_mat & symmSparseMatrix) {
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

sp_mat AugTree::ComputeHpred(const mat & spCoords, const vec & time, const mat & covariateMatrix) {

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

  return updateHmatrixPred() ;
}

vec AugTree::ComputeEvar(const sp_mat & HmatPred, const int batchSize) {

  double errorVar = std::pow(m_errorSD, 2) ;
  vec EvarValues(HmatPred.n_rows, fill::zeros) ;
  int obsIndex = 0 ;

  while (obsIndex < HmatPred.n_rows) {

    int newObsIndex = std::min(obsIndex + batchSize - 1, int(HmatPred.n_rows) - 1) ;
    sp_mat bVecTrans = HmatPred.block(obsIndex, 0, newObsIndex - obsIndex + 1, HmatPred.cols()).transpose() ;
    sp_mat meanValue = bVecTrans % m_FullCondPrecisionChol.solve(bVecTrans) ;
    EvarValues.subvec(obsIndex, newObsIndex) = (sum(meanValue,0)).transpose() + errorVar ;

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

  bool test = (m_spacetimeScaling == scalePara) && (m_MRAcovParasSpace == MRAcovParasSpace) && (m_MRAcovParasTime == MRAcovParasTime) ;

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
