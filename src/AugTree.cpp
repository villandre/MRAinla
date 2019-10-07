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
using namespace std ;

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

AugTree::AugTree(uint & M, Array2d & lonRange, Array2d & latRange, Array2d & timeRange, vec & observations, ArrayXXd & obsSp, ArrayXd & obsTime, ArrayXXd & predCovariates, ArrayXXd & predSp, ArrayXd & predTime, uint & minObsForTimeSplit, unsigned long int & seed, ArrayXXd & covariates, const bool splitTime, const unsigned int numKnotsRes0, const unsigned int J)
  : m_M(M)
{
  m_GammaParasSet = false ;
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

  createHmatrixPred() ;
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

// In this version, we first split the longitude, then the latitude, and finally time.
// We make sure that splits don't result in empty regions
void AugTree::createLevels(TreeNode * parent, const uint & numObsForTimeSplit, const bool splitTime) {
  ArrayXi obsForMedian = parent->GetObsInNode() ;
  ArrayXi childMembership = ArrayXi::Zero(obsForMedian.size()) ;
  std::vector<dimensions> childDimensions ;
  childDimensions.push_back(parent->GetDimensions()) ;

  if (obsForMedian.size() <= 1) {
    throw Rcpp::exception("Cannot split empty region or region with only one observation.\n") ;
  }

  // Handling longitude

  ArrayXi elementsInChild = ArrayXi::LinSpaced(obsForMedian.size(), 0, obsForMedian.size()-1) ;
  ArrayXd column = m_dataset.spatialCoords.col(0) ;
  ArrayXd elementsForMedian = elem(column, obsForMedian) ;
  double colMedian = median(elementsForMedian) ;
  ArrayXd updatedLongitude(2) ;
  updatedLongitude(0) = childDimensions.at(0).longitude(0) ;
  updatedLongitude(1) = colMedian ;
  ArrayXd newChildLongitude(2) ;
  newChildLongitude(0) = colMedian ;
  newChildLongitude(1) = childDimensions.at(0).longitude(1) ;
  dimensions newDimensions = childDimensions.at(0) ;
  newDimensions.longitude = newChildLongitude ;
  childDimensions.push_back(newDimensions) ;
  childDimensions.at(0).longitude = updatedLongitude ;
  ArrayXi greaterElements = find(elementsForMedian > colMedian) ; // In deriveObsInNode, the checks are <=. It follows that observations on the right and upper boundaries of a zone are included.
  ArrayXi updateIndices = elem(elementsInChild, greaterElements) ;
  for (uint innerIndex = 0; innerIndex < updateIndices.size(); innerIndex++) {
    childMembership(updateIndices(innerIndex)) = 1 ;
  }

  uint newChildIndex = 2 ;

  // Handling latitude

  for (int j = 0 ; j < 2; j++) {
    ArrayXd column = m_dataset.spatialCoords.col(1) ;
    ArrayXd subColumn = elem(column, obsForMedian) ;
    ArrayXi elementsInChild = find(childMembership == j) ;
    double colMedian = median(elem(subColumn, elementsInChild)) ;
    ArrayXd updatedLatitude(2) ;
    updatedLatitude(0) = childDimensions.at(j).latitude(0);
    updatedLatitude(1) = colMedian ;
    ArrayXd newChildLatitude(2) ;
    newChildLatitude(0) = colMedian ;
    newChildLatitude(1) = childDimensions.at(j).latitude(1) ;
    dimensions newDimensions = childDimensions.at(j) ;
    newDimensions.latitude = newChildLatitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(j).latitude = updatedLatitude ;
    ArrayXi greaterElements = find(elem(subColumn, elementsInChild).array() > colMedian) ;
    ArrayXi updateIndices = elem(elementsInChild, greaterElements) ;
    for (uint innerIndex = 0; innerIndex < updateIndices.size(); innerIndex++) {
      childMembership(updateIndices(innerIndex)) = newChildIndex ;
    }
    newChildIndex += 1 ;
  }
  if (splitTime) {
    // Handling time
    for (int j = 0 ; j < newChildIndex; j++) {
      ArrayXd time = elem(m_dataset.timeCoords, obsForMedian) ;
      ArrayXi elementsInChild = find(childMembership == j) ;
      ArrayXd elementsForMedian = elem(time, elementsInChild) ;
      bool onlyOneTimePoint = true;
      uint innerIndex = 1 ;
      while(onlyOneTimePoint) {
        onlyOneTimePoint = (elementsForMedian(innerIndex) == elementsForMedian(0)) ;
        innerIndex += 1 ;
        if (innerIndex == elementsForMedian.size()) break ;
      }

      if (!onlyOneTimePoint) {
        double colMedian = median(elementsForMedian) ;
        ArrayXd updatedTime(2) ;
        updatedTime(0)  = childDimensions.at(j).time(0) ;
        updatedTime(1) = colMedian ;
        ArrayXd newChildTime(2) ;
        newChildTime(0) = colMedian ;
        newChildTime(1) = childDimensions.at(j).time(1) ;
        dimensions newDimensions = childDimensions.at(j) ;
        newDimensions.time = newChildTime ;
        childDimensions.push_back(newDimensions) ;
        childDimensions.at(j).time = updatedTime ;
        ArrayXi greaterElements = find(elem(time, elementsInChild).array() > colMedian) ;
        ArrayXi updateIndices = elem(elementsInChild, greaterElements) ;
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

void AugTree::CreateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  mat FEinvMatrix = pow(m_fixedEffSD, -2) * mat::Identity(m_fixedEffParameters.size(), m_fixedEffParameters.size()) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;

  // This is because we want the Q matrix to have a block-diagonal component in the lower-right corner,
  // which prompted us to make H* = [X,F] rather than H* = [F, X], like in Eidsvik.

  m_SigmaFEandEtaInv = createBlockMatrix(FEinvAndKinvMatrixList) ;
}

void AugTree::UpdateSigmaBetaEtaInvMat() {
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

  for (int k = 0; k < m_FEinvAndKinvMatrices.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_FEinvAndKinvMatrices, k); it; ++it)
    {
      it.valueRef() = concatenatedValues(loopIndex) ;
      loopIndex += 1 ;
    }
  }
}

void AugTree::createHmatrix() {

  cout << "Creating the data mapping matrix... \n" ;
  int numObs = m_dataset.spatialCoords.rows() ;
  ArrayXXi HmatPos(0, 2) ;

  ArrayXi FmatObsOrder = ArrayXi::Zero(numObs) ;
  int rowIndex = 0 ;
  int colIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
  int numTips = tipNodes.size() ;

  std::vector<ArrayXi> ancestorIdsVec(numTips) ;

  for (uint i = 0 ; i < tipNodes.size(); i++) {
    ArrayXi idVec = tipNodes.at(i)->GetAncestorIds() ; // Last element is tip node.
    ancestorIdsVec.at(i) = idVec ;
  }

  std::sort(tipNodes.begin(), tipNodes.end(), [] (TreeNode * first, TreeNode * second) {
    ArrayXi firstAncestorIds = first->GetAncestorIds() ;
    ArrayXi secondAncestorIds = second->GetAncestorIds() ;
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
    ArrayXi observationIndices = nodeToProcess->GetObsInNode() ;

    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.segment(rowIndex, observationIndices.size()) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      ArrayXi brickAncestors(brickList.size()) ;
      for (uint i = 0 ; i < brickAncestors.size() ; i++) {
        brickAncestors(i) = brickList.at(i)->GetNodeId() ;
      }
      uint index = 0 ;
      for (uint nodeIndex = 0; nodeIndex < m_vertexVector.size() ; nodeIndex++) {
        int nodePos = GetNodePos(nodeIndex) ;
        if (index < brickAncestors.size()) {
          if (nodeIndex == brickAncestors(index)) {
            ArrayXi rowIndices = rep(uvec::LinSpaced(nodeToProcess->GetB(index).rows(), 0, nodeToProcess->GetB(index).rows() - 1).array(),
                                  nodeToProcess->GetB(index).cols()) + rowIndex  ;
            ArrayXi colIndices = rep_each(uvec::LinSpaced(nodeToProcess->GetB(index).cols(), 0, nodeToProcess->GetB(index).cols() - 1).array(),
                                       nodeToProcess->GetB(index).rows()) + colIndex ;

            ArrayXXi positions = join_rows(rowIndices, colIndices) ;
            HmatPos = join_cols(HmatPos, positions) ; // Slow, but this is not done many times.

            index += 1 ;
          }
        }
        colIndex += m_vertexVector.at(nodePos)->GetNumKnots() ;
      }
      colIndex = 0 ;
      rowIndex += observationIndices.size() ; // The B matrices should have as many rows as observations in the node...
    }
  }

  HmatPos.col(1) += m_dataset.covariateValues.cols() + 1 ;
  m_obsOrderForFmat = FmatObsOrder ;

  mat intercept = mat::Ones(numObs, 1) ;
  mat incrementedCovar = join_rows(intercept.array(), m_dataset.covariateValues.array()) ;
  MatrixXd incrementedCovarReshuffled = rows(incrementedCovar.array(), m_obsOrderForFmat) ;

  // The idea behind this is to associate a memory location in the W matrices and a cell
  // in the H matrix. The memory location we get from calling .data() on a given W
  // matrix changes when the W values are updated, even when W is allocated on the
  // heap, i.e. with 'new', but the memory location of the vector element is constant
  // once the vector is fully populated. To get access to the value, all we need to do
  // is invoke *(pointer->data() + offset).
  // The odd nesting of the loops is due to the fact that m_Hmat is row-major.
  // I made m_Hmat row-major to make the updating process, as all values on any
  // given row are associated with the same tip.
  // Don't mind the four nested loops: all it does is traverse all the elements once.
  for (auto & tipNode : tipNodes) {
    uint numObs = tipNode->GetObsInNode().size() ;
    for (uint rowIndex = 0; rowIndex < numObs ; rowIndex++) {
      for (uint depth = 0 ; depth <= m_M ; depth++) {
        uint numKnotsAtDepth = tipNode->GetWlist().at(depth).cols() ;
        for (uint colIndex = 0; colIndex < numKnotsAtDepth; colIndex++) {
          pointerOffset elementToAdd = pointerOffset(&(tipNode->GetWlist().at(depth)), colIndex * numObs + rowIndex) ;
          m_pointerOffsetForHmat.push_back(elementToAdd) ;
        }
      }
    }
  }

  ArrayXXi covPos = join_rows(rep(uvec::LinSpaced(incrementedCovarReshuffled.rows(), 0, incrementedCovarReshuffled.rows() - 1).array(), incrementedCovarReshuffled.cols()),
                          rep_each(uvec::LinSpaced(incrementedCovarReshuffled.cols(), 0, incrementedCovarReshuffled.cols() - 1).array(), incrementedCovarReshuffled.rows())) ;
  HmatPos = join_cols(covPos, HmatPos).eval() ; // Could this cause aliasing?

  ArrayXd concatenatedValues = ArrayXd::Zero(HmatPos.rows()) ;
  concatenatedValues.segment(0, incrementedCovarReshuffled.size()) = incrementedCovarReshuffled ;
  uint secondIndex = incrementedCovar.size() ;
  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetObsInNode().size() > 0) {
      for (auto & Bmat : tipNode->GetWlist()) {
        concatenatedValues.segment(secondIndex, Bmat.size()) = Bmat ;
        secondIndex += Bmat.size() ;
      }
    }
  }
  m_Hmat.resize(incrementedCovarReshuffled.rows(), m_numKnots + incrementedCovar.cols()) ;
  std::vector<Triplet> tripletList;

  int offset = 0 ;

  for (uint i = 0; i < concatenatedValues.size(); i++) {
    tripletList.push_back(Triplet(HmatPos(i, 0), HmatPos(i, 1), concatenatedValues(i))) ;
  }

  m_Hmat.setFromTriplets(tripletList.begin(), tripletList.end()) ;
  std::cout << "Done! \n" ;
}

void AugTree::updateHmatrix() {
  int loopIndex = 0 ;
  for (int k=0; k<m_Hmat.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_Hmat, k); it; ++it)
    {
      if (it.col() >= m_fixedEffParameters.size()) { // Covariates never change, only the elements in the F matrix change.
        it.valueRef() = *(m_pointerOffsetForHmat.at(loopIndex).matrixLocation->data() + m_pointerOffsetForHmat.at(loopIndex).offset) ;
        loopIndex += 1 ;
      }
    }
  }
}

void AugTree::updateHmatrixPred() {
  mat intercept = ArrayXd::Ones(m_predictData.covariateValues.rows(), 1) ;
  ArrayXXd incrementedCovarPredict = join_rows(intercept.array(), m_predictData.covariateValues.array()) ;
  ArrayXXd incrementedCovarPredictReshuffled = rows(incrementedCovarPredict, m_obsOrderForHpredMat) ;
  // Will this correctly preserve the order? It's supposed to...
  ArrayXd concatenatedValues(Map<ArrayXd>(incrementedCovarPredictReshuffled.data(), incrementedCovarPredictReshuffled.cols() * incrementedCovarPredictReshuffled.rows()))  ;
  int totalSize = concatenatedValues.size() ;
  int secondIndex = totalSize ;
  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetPredIndices().size() > 0) {
      for (auto & Umat : tipNode->GetUmatList()) {
        totalSize += Umat.size();
      }
    }
  }
  concatenatedValues.conservativeResize(totalSize, 1) ;

  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetPredIndices().size() > 0) {
      for (auto & Umat : tipNode->GetUmatList()) {
        ArrayXd B(Map<ArrayXd>(Umat.data(), Umat.cols() * Umat.rows()));
        concatenatedValues.segment(secondIndex, B.size()) = B ;
        secondIndex += B.size() ;
      }
    }
  }

  int loopIndex = 0 ;
  for (int k = 0; k < m_HmatPred.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_HmatPred, k); it; ++it)
    {
      it.valueRef() = concatenatedValues(loopIndex) ;
      loopIndex += 1 ;
    }
  }
}

void AugTree::createHmatrixPred() {
  std::cout << "Creating the prediction mapping matrix H... \n" ;
  uint numObs = m_predictData.spatialCoords.rows() ;
  ArrayXXi HmatPredPos(0, 2) ;

  ArrayXi FmatObsOrder = ArrayXi::Zero(numObs) ;
  uint rowIndex = 0 ;
  uint colIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
  int numTips = tipNodes.size() ;

  std::vector<ArrayXi> ancestorIdsVec(numTips) ;

  for (uint i = 0 ; i < tipNodes.size(); i++) {
    ArrayXi idVec = tipNodes.at(i)->GetAncestorIds() ; // Last element is tip node.
    ancestorIdsVec.at(i) = idVec ;
  }

  std::sort(tipNodes.begin(), tipNodes.end(), [] (TreeNode * first, TreeNode * second) {
    ArrayXi firstAncestorIds = first->GetAncestorIds() ;
    ArrayXi secondAncestorIds = second->GetAncestorIds() ;
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

  m_obsOrderForHpredMat = ArrayXi::Zero(m_predictData.timeCoords.size()) ;
  uint elementIndex = 0 ;
  for (auto & i : tipNodes) {
    uint numPredObsInNode = i->GetPredIndices().size() ;
    if (numPredObsInNode > 0) {
      m_obsOrderForHpredMat.segment(elementIndex, i->GetPredIndices().size()) = i->GetPredIndices() ;
      elementIndex += i->GetPredIndices().size() ;
    }
  }

  for (auto & nodeToProcess : tipNodes) {
    ArrayXi observationIndices = nodeToProcess->GetPredIndices() ;
    std::vector<TreeNode *> ancestorsList = nodeToProcess->getAncestors() ;

    bool test =  observationIndices.size() > 0 ;

    if (test) {
      FmatObsOrder.segment(rowIndex, observationIndices.size()) = observationIndices ;

      std::vector<TreeNode *> brickList = nodeToProcess->getAncestors() ;

      ArrayXi brickAncestors(brickList.size()) ;
      for (uint i = 0 ; i < brickAncestors.size() ; i++) {
        brickAncestors(i) = brickList.at(i)->GetNodeId() ;
      }
      uint index = 0 ;

      for (uint nodeIndex = 0; nodeIndex < m_vertexVector.size() ; nodeIndex++) {
        int nodePos = GetNodePos(nodeIndex) ;
        if (index < brickAncestors.size()) {
          if (nodeIndex == brickAncestors(index)) {
            ArrayXi rowIndices = rep(uvec::LinSpaced(observationIndices.size(), 0, observationIndices.size() - 1).array(),
                                  ancestorsList.at(index)->GetNumKnots()) + rowIndex  ;
            ArrayXi colIndices = rep_each(uvec::LinSpaced(ancestorsList.at(index)->GetNumKnots(), 0, ancestorsList.at(index)->GetNumKnots() - 1).array(),
                                       observationIndices.size()) + colIndex  ;
            ArrayXXi positions = join_rows(rowIndices, colIndices) ;
            HmatPredPos = join_cols(HmatPredPos, positions) ; // Slow, but this is not done many times.

            index += 1 ;
          }
        }
        colIndex += m_vertexVector.at(nodePos)->GetNumKnots() ;
      }
      colIndex = 0 ;
      rowIndex += observationIndices.size() ;
    }
  }

  mat intercept = mat::Ones(numObs, 1) ;
  mat transIncrementedCovar = (join_rows(intercept.array(), m_predictData.covariateValues.array())).transpose() ;
  mat incrementedCovarPredictReshuffled = cols(transIncrementedCovar.array(), FmatObsOrder).transpose() ;

  HmatPredPos.col(1) += m_predictData.covariateValues.cols() + 1 ;
  ArrayXXi covPos = join_rows(rep(uvec::LinSpaced(m_predictData.covariateValues.rows(), 0, m_predictData.covariateValues.rows() - 1).array(),
                                  m_predictData.covariateValues.cols() + 1),
                          rep_each(uvec::LinSpaced(m_predictData.covariateValues.cols() + 1, 0, m_predictData.covariateValues.cols()).array(),
                                   m_predictData.covariateValues.rows())) ;
  HmatPredPos = join_cols(covPos, HmatPredPos).eval() ; // Could this cause aliasing?

  ArrayXd concatenatedValues(Map<ArrayXd>(incrementedCovarPredictReshuffled.data(), incrementedCovarPredictReshuffled.cols() * incrementedCovarPredictReshuffled.rows()))  ;
  int totalSize = concatenatedValues.size() ;
  int secondIndex = totalSize ;
  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetPredIndices().size() > 0) {
      for (auto & Umat : tipNode->GetUmatList()) {
        totalSize += Umat.size();
      }
    }
  }
  concatenatedValues.conservativeResize(totalSize, 1) ;

  for (auto & tipNode : GetTipNodes()) {
    if (tipNode->GetPredIndices().size() > 0) {
      for (auto & Umat : tipNode->GetUmatList()) {
        ArrayXd B(Map<ArrayXd>(Umat.data(), Umat.cols() * Umat.rows()));
        concatenatedValues.segment(secondIndex, B.size()) = B ;
        secondIndex += B.size() ;
      }
    }
  }
  std::vector<Triplet> tripletList;

  int offset = 0 ;
  m_HmatPred.resize(m_predictData.timeCoords.size(), m_numKnots + incrementedCovarPredictReshuffled.cols()) ;

  for (uint i = 0; i < concatenatedValues.size(); i++) {
    tripletList.push_back(Triplet(HmatPredPos(i, 0), HmatPredPos(i, 1), concatenatedValues(i))) ;
  }
  m_HmatPred.setFromTriplets(tripletList.begin(), tripletList.end()) ;
  std::cout << "Done! \n" ;
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
  fflush(stdout);
  int n = m_dataset.responseValues.size() ;

  if (m_SigmaFEandEtaInv.rows() == 0) {
    CreateSigmaBetaEtaInvMat() ;
  } else {
    UpdateSigmaBetaEtaInvMat() ;
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> blockChol(m_SigmaFEandEtaInv) ;

  double logDetSigmaKFEinv = blockChol.vectorD().array().log().sum() ;

  if (m_recomputeMRAlogLik) {
    if (m_Hmat.size() == 0) {
      createHmatrix() ;
    } else {
      updateHmatrix() ;
    }
  }

  sp_mat secondTerm ;
  secondTerm = std::pow(m_errorSD, -2) * m_Hmat.transpose() * m_Hmat ;
  vec responsesReshuffled = elem(m_dataset.responseValues.array(), m_obsOrderForFmat) ;

  mat scaledResponse ;
  scaledResponse.noalias() = std::pow(m_errorSD, -2) * responsesReshuffled.transpose() * m_Hmat ;

  m_FullCondPrecisionChol.compute(m_SigmaFEandEtaInv + secondTerm) ;

  // m_FullCondSDs = sqrt(m_FullCondPrecisionChol.diag()) ;
  // cout << "Done... \n" ;

  vec updatedMean = m_FullCondPrecisionChol.solve(scaledResponse.transpose()) ;
  if(m_FullCondPrecisionChol.info()!=Success) {
    // solving failed
    std::cout<< "Solving failed!!!! \n" ;
    throw Rcpp::exception("Leave now... \n") ;
  }

  m_Vstar = updatedMean ; // Assuming there will be an implicit conversion to vec type.

  m_FullCondMean = m_Vstar ;

  vec fixedEffMeans = m_Vstar.head(m_fixedEffParameters.size()) ;
  SetFixedEffParameters(fixedEffMeans) ;

  double logDetQmat = m_FullCondPrecisionChol.vectorD().array().log().sum() ; // In LDLT, the L matrix has a diagonal with 1s, meaning that its determinant is 1. It follows that the determinant of the original matrix is simply the product of the elements in the D matrix.

  m_logFullCond = 0.5 * logDetQmat ; // Since we arbitrarily evaluate always at the full-conditional mean, the exponential part of the distribution reduces to 0.

  // Computing p(v* | Psi)
  mat logLikTerm ;
  logLikTerm.noalias() = -0.5 * m_Vstar.transpose() * m_SigmaFEandEtaInv.selfadjointView<Lower>() * m_Vstar ;

  m_logCondDist = logLikTerm(0,0) + 0.5 * logDetSigmaKFEinv ;

  // Computing p(y | v*, Psi)
  double errorLogDet = -2 * n * log(m_errorSD) ;
  vec recenteredY ;
  recenteredY.noalias() = responsesReshuffled - m_Hmat * m_Vstar ;

  vec globalLogLikExp ;
  globalLogLikExp.noalias() = -0.5 * std::pow(m_errorSD, -2) * recenteredY.transpose() * recenteredY ;
  m_globalLogLik = 0.5 * errorLogDet + globalLogLikExp(0) ;
}

void AugTree::ComputeLogJointPsiMarginal() {

  ComputeLogPriors() ;

  if (m_recomputeMRAlogLik) {
    m_MRAcovParasSpace.print("Space parameters:") ;
    m_MRAcovParasTime.print("Time parameters:") ;
    computeWmats() ; // This will produce the K matrices required.
  }

  ComputeLogFCandLogCDandDataLL() ;

  printf("Observations log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
   m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  m_logJointPsiMarginal = m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond ;
  printf("Joint value: %.4e \n \n", m_logJointPsiMarginal) ;
}

void AugTree::ComputeHpred(const mat & spCoords, const vec & time, const mat & covariateMatrix) {

  // SetPredictData(spCoords, time, covariateMatrix) ;
  std::vector<TreeNode *> tipNodes = GetTipNodes() ;

  for (auto & i : tipNodes) {
    // i->SetPredictLocations(m_predictData) ;
    ArrayXi predictionsInLeaf = i->GetPredIndices() ;
    if (predictionsInLeaf.size() > 0) {
      i->computeUpred(m_MRAcovParasSpace, m_MRAcovParasTime, m_spacetimeScaling, m_predictData, m_matern, m_spaceNuggetSD, m_timeNuggetSD) ;
    }
  }
  if (m_HmatPred.rows() == 0) {
    createHmatrixPred() ;
  } else {
    updateHmatrixPred() ;
  }
}

vec AugTree::ComputeEvar(const int batchSize) {

  double errorVar = std::pow(m_errorSD, 2) ;
  vec EvarValues = vec::Zero(m_HmatPred.rows()) ;
  int obsIndex = 0 ;

  while (obsIndex < m_HmatPred.rows()) {

    int newObsIndex = std::min(obsIndex + batchSize - 1, int(m_HmatPred.rows()) - 1) ;
    mat bVecTrans = m_HmatPred.block(obsIndex, 0, newObsIndex - obsIndex + 1, m_HmatPred.cols()).transpose() ;
    mat meanValue = m_FullCondPrecisionChol.solve(bVecTrans) ;
    EvarValues.segment(obsIndex, newObsIndex - obsIndex + 1) = meanValue.colwise().sum().transpose() + errorVar * mat::Ones(newObsIndex - obsIndex + 1, 1) ;

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
