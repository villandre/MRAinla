// #ifdef _OPENMP
#include <omp.h>
// #endif

#include <math.h>

#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"
// #include "gperftools/profiler.h"

using namespace Rcpp ;
using namespace MRAinla ;
using namespace Eigen ;
using namespace std ;

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

AugTree::AugTree(const uint & Mlon,
                 const uint & Mlat,
                 const uint & Mtime,
                 const Array2d & lonRange,
                 const Array2d & latRange,
                 const Array2d & timeRange,
                 const vec & observations,
                 const ArrayXXd & obsSp,
                 const ArrayXd & obsTime,
                 const ArrayXXd & predCovariates,
                 const ArrayXXd & predSp,
                 const ArrayXd & predTime,
                 const uint & minObsForTimeSplit,
                 const unsigned long int & seed,
                 const ArrayXXd & covariates,
                 const bool & splitTime,
                 const unsigned int & numKnotsRes0,
                 const double & J,
                 const string & distMethod,
                 const Rcpp::List & STparsHyperpars,
                 const Rcpp::NumericVector & fixedEffParsHyperpars,
                 const Rcpp::NumericVector & errorParsHyperpars,
                 const Rcpp::NumericVector & FEmuVec,
                 const double & nuggetSD,
                 const bool & normalHyperprior)
  : m_distMethod(distMethod), m_nuggetSD(nuggetSD),
    m_Mlon(Mlon), m_Mlat(Mlat), m_Mtime(Mtime),
    m_normalHyperprior(normalHyperprior)
{
  m_M = Mlon + Mlat + Mtime ;
  m_dataset = inputdata(observations, obsSp, obsTime, covariates) ;
  m_mapDimensions = dimensions(lonRange, latRange, timeRange) ;
  std::random_device rd ;
  std::mt19937_64 generator(rd()) ;
  generator.seed(seed) ;
  m_randomNumGenerator = generator ;
  SetPredictData(predSp, predTime, predCovariates) ;
  m_assignedPredToKnot = Array<bool, Dynamic, 1>(m_predictData.timeCoords.size()) ;
  m_assignedPredToKnot.segment(0, m_assignedPredToKnot.size()) = false ;

  m_fixedEffParameters = Eigen::VectorXd::Zero(m_dataset.covariateValues.cols() + 1);
  Rcout << "Building the grid..." << std::endl ;
  BuildTree(minObsForTimeSplit, splitTime, numKnotsRes0, J) ;
  Rcout << "Grid built!" << std::endl ;

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

  SetParsHyperpars(STparsHyperpars, fixedEffParsHyperpars, errorParsHyperpars, normalHyperprior) ;

  m_FEmu = Rcpp::as<vec>(FEmuVec) ;
  // for (auto & i : m_vertexVector) {
  //   Rprintf("This is node %i. I contain %i knots and %i observations. \n", i->GetNodeId(), i->GetNumKnots(), i->GetNumObs()) ;
  // }
}

void AugTree::BuildTree(const uint & minObsForTimeSplit, const bool splitTime, const unsigned int numKnots0, double J)
{
  m_vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(m_mapDimensions, m_dataset) ;

  m_vertexVector.push_back(topNode) ;
  ArrayXi numSplits(3) ;
  numSplits << m_Mlon, m_Mlat, m_Mtime;
  createLevels(topNode, "longitude", numSplits) ;
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

// We make sure that splits don't result in empty regions
// numObsForTimeSplit still does nothing

void AugTree::createLevels(TreeNode * parent, std::string splitWhat, ArrayXi numSplitsLeft) {
  ArrayXi obsForMedian = parent->GetObsInNode() ;
  ArrayXi childMembership = ArrayXi::Zero(obsForMedian.size()) ;
  std::vector<dimensions> childDimensions ;
  childDimensions.push_back(parent->GetDimensions()) ;

  if (obsForMedian.size() <= 1) {
    throw Rcpp::exception("Cannot split empty region or region with only one observation.\n") ;
  }

  // Handling longitude
  if (splitWhat == "longitude") {
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
  } else if (splitWhat == "latitude") {
    ArrayXi elementsInChild = ArrayXi::LinSpaced(obsForMedian.size(), 0, obsForMedian.size()-1) ;
    ArrayXd column = m_dataset.spatialCoords.col(1) ;
    ArrayXd elementsForMedian = elem(column, obsForMedian) ;
    double colMedian = median(elementsForMedian) ;
    ArrayXd updatedLatitude(2) ;
    updatedLatitude(0) = childDimensions.at(0).latitude(0) ;
    updatedLatitude(1) = colMedian ;
    ArrayXd newChildLatitude(2) ;
    newChildLatitude(0) = colMedian ;
    newChildLatitude(1) = childDimensions.at(0).latitude(1) ;
    dimensions newDimensions = childDimensions.at(0) ;
    newDimensions.latitude = newChildLatitude ;
    childDimensions.push_back(newDimensions) ;
    childDimensions.at(0).latitude = updatedLatitude ;
  } else if (splitWhat == "time") {

    bool onlyOneTimePoint = true;
    uint innerIndex = 1 ;
    ArrayXd elementsForMedian = elem(m_dataset.timeCoords, obsForMedian) ;
    while(onlyOneTimePoint) {
      onlyOneTimePoint = (elementsForMedian(innerIndex) - elementsForMedian(0)) < 3e-3  ; // This is to account for the jittering
      innerIndex += 1 ;
      if (innerIndex == elementsForMedian.size()) break ;
    }
    if (!onlyOneTimePoint) {
      double colMedian = median(elementsForMedian) ;
      ArrayXd updatedTime(2) ;
      updatedTime(0)  = childDimensions.at(0).time(0) ;
      updatedTime(1) = colMedian ;
      ArrayXd newChildTime(2) ;
      newChildTime(0) = colMedian ;
      newChildTime(1) = childDimensions.at(0).time(1) ;
      dimensions newDimensions = childDimensions.at(0) ;
      newDimensions.time = newChildTime ;
      childDimensions.push_back(newDimensions) ;
      childDimensions.at(0).time = updatedTime ;
    }
  }

  int incrementedDepth = parent->GetDepth() + 1 ;
  for (auto & i : childDimensions) {
    if (parent->GetDepth() < m_M - 1) {
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
  // Domain is split with respect to longitude, then latitude, then time.
  if (incrementedDepth < m_M) {
    if (splitWhat == "longitude") {
      if (numSplitsLeft(1) > 0) {
        splitWhat = "latitude" ;
        numSplitsLeft(1) -= 1 ;
      } else if (numSplitsLeft(2) > 0){
        splitWhat = "time" ;
        numSplitsLeft(2) -= 1 ;
      } else {
        numSplitsLeft(0) -= 1 ;
      }
    } else if (splitWhat == "latitude") {
      if (numSplitsLeft(2) > 0){
        splitWhat = "time" ;
        numSplitsLeft(2) -= 1 ;
      } else if (numSplitsLeft(0) > 0) {
        splitWhat = "longitude" ;
        numSplitsLeft(0) -= 1 ;
      } else {
        numSplitsLeft(1) -= 1 ;
      }
    } else {
      if (numSplitsLeft(0) > 0) {
        splitWhat = "longitude" ;
        numSplitsLeft(0) -= 1 ;
      } else if (numSplitsLeft(1) > 0) {
        splitWhat = "latitude" ;
        numSplitsLeft(1) -= 1 ;
      } else {
        numSplitsLeft(2) -= 1 ;
      }
    }
    for (auto && i : parent->GetChildren()) {
      createLevels(i, splitWhat, numSplitsLeft) ;
    }
  }
}

void AugTree::generateKnots(TreeNode * node, const unsigned int numKnotsRes0, double J) {

  int numNodesAtLevel = GetLevelNodes(node->GetDepth()).size() ;
  int numKnotsToGen = std::max(uint(std::ceil((numKnotsRes0 * pow(J, node->GetDepth()))/numNodesAtLevel)), uint(2)) ;
  // node->genRandomKnots(m_dataset, numKnotsToGen, m_randomNumGenerator) ;
  if (node->GetDepth() < m_M) {
    node->genKnotsOnCube(m_predictData, numKnotsToGen, m_randomNumGenerator, m_assignedPredToKnot) ;
  } else {
    node->genKnotsOnCube(m_dataset, numKnotsToGen, m_randomNumGenerator, m_assignedPredToKnot) ;
  }
  if (node->GetChildren().at(0) != NULL) {
    for (auto &i : node->GetChildren()) {
      generateKnots(i, numKnotsRes0, J) ;
    }
  }
}

void AugTree::computeWmats() {
  m_vertexVector.at(0)->ComputeWmat(m_MaternParsSpace, m_MaternParsTime, m_spacetimeScaling, m_nuggetSD, m_distMethod) ;

  for (uint level = 1; level <= m_M; level++) {
    std::vector<TreeNode *> levelNodes = GetLevelNodes(level) ;

    // Trying openmp. We need to have a standard looping structure.
    // pragma omp parallel for

    for (std::vector<TreeNode *>::iterator it = levelNodes.begin(); it < levelNodes.end(); it++)
    {
      (*it)->ComputeWmat(m_MaternParsSpace, m_MaternParsTime, m_spacetimeScaling, m_nuggetSD, m_distMethod) ;
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

void AugTree::CreateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  // We have to use that loop instead of creating an identity matrix, else, it will be assumed in the updating function
  // that the zeros in the identity matrix are not ALWAYS 0 in the sparse matrix, and the iterating across
  // non-zero elements will cover those zeros too!
  mat FEinvMatrix(1,1) ;
  FEinvMatrix(0,0) = pow(m_fixedEffSD, -2) ;
  for (uint i = 0; i < m_fixedEffParameters.size(); i++) {
    FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;
  }

  m_SigmaFEandEtaInv = createBlockMatrix(FEinvAndKinvMatrixList) ;
}

void AugTree::UpdateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList = getKmatricesInversePointers() ;
  // We have to use that loop instead of creating an identity matrix, else, it will be assumed in the updating function
  // that the zeros in the identity matrix are not ALWAYS 0 in the sparse matrix, and the iterating across
  // non-zero elements will cover those zeros too!
  mat FEinvMatrix(1,1) ;
  FEinvMatrix(0,0) = pow(m_fixedEffSD, -2) ;
  for (uint i = 0; i < m_fixedEffParameters.size(); i++) {
    FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &FEinvMatrix) ;
  }

  sp_mat::InnerIterator it(m_SigmaFEandEtaInv, 0) ;
  int k = 0 ;

  for (auto & matrixPointer : FEinvAndKinvMatrixList) {
    int colIndex = 0 ;
    for (uint elementIndex = 0 ; elementIndex < matrixPointer->size(); elementIndex++) {
      it.valueRef() = (*matrixPointer)(elementIndex) ;
      colIndex += 1 ;
      if (colIndex == matrixPointer->cols()) {
        it = sp_mat::InnerIterator(m_SigmaFEandEtaInv, k + 1) ;
        k += 1 ;
        colIndex = 0 ;
      } else {
        it = ++it ;
      }
    }
  }
}

void AugTree::createHmatrix() {
  int numObs = m_dataset.spatialCoords.rows() ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
  int numTips = tipNodes.size() ;

  ArrayXi FmatObsOrder = ArrayXi::Zero(numObs) ;
  int rowIndex = 0 ;
  std::vector<ArrayXi> ancestorIdsVec(numTips) ;

  for (uint i = 0 ; i < tipNodes.size(); i++) {
    ArrayXi idVec = tipNodes.at(i)->GetAncestorIds() ; // Last element is tip node.
    ancestorIdsVec.at(i) = idVec ;
  }

  std::vector<Triplet> tripletList ;

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

  ArrayXi colIndexAtEachRes = ArrayXi::Zero(m_M + 1);

  for (auto & node : m_vertexVector) {
    if (node->GetDepth() < m_M) {
      for (uint i = node->GetDepth() + 1; i <= m_M; i++) {
        colIndexAtEachRes(i) += node->GetNumKnots() ;
      }
    }
  } // Convoluted way of getting cumsum(c(0, Num knots at resolutions 0 to M-1))

  colIndexAtEachRes += (m_dataset.covariateValues.cols() + 1) ; // Covariate columns go before, the +1 is for the intercept.
  std::vector<TreeNode *> previousBrickAncestors = tipNodes.at(0)->getAncestors() ;

  for (auto & nodeToProcess : tipNodes) {

    if (nodeToProcess->GetObsInNode().size() == 0 ) {
      Rcpp::stop("Why is there an empty node? This is not supposed to happen. \n\n") ;
    } // Should not normally occur, as the tree-building scheme forbids empty nodes.

    // The idea behind this section of the code is that the column index for producing the section
    // of the H matrix for a tip node at any given depth i should only change when a different ancestor for the
    // tip node is reached. The hierarchical structure of the tree explains this.
    for (uint i = 1 ; i <= m_M; i++) { // The root node is shared by all tips, hence the 1.
      if (previousBrickAncestors.at(i) != nodeToProcess->getAncestors().at(i)) {
        colIndexAtEachRes(i) += previousBrickAncestors.at(i)->GetNumKnots() ;
      }
    }
    previousBrickAncestors = nodeToProcess->getAncestors() ;

    FmatObsOrder.segment(rowIndex, nodeToProcess->GetObsInNode().size()) = nodeToProcess->GetObsInNode() ;

    for (uint depthIndex = 0; depthIndex <= m_M; depthIndex++) {
      ArrayXi rowIndices = rep(uvec::LinSpaced(nodeToProcess->GetB(depthIndex).rows(), 0, nodeToProcess->GetB(depthIndex).rows() - 1).array(),
          nodeToProcess->GetB(depthIndex).cols()) + rowIndex  ;
      ArrayXi colIndices = rep_each(uvec::LinSpaced(nodeToProcess->GetB(depthIndex).cols(), 0, nodeToProcess->GetB(depthIndex).cols() - 1).array(),
          nodeToProcess->GetB(depthIndex).rows()) + colIndexAtEachRes(depthIndex) ;

      for (uint i = 0; i < rowIndices.size(); i++) {
        tripletList.push_back(Triplet(rowIndices(i), colIndices(i), nodeToProcess->GetB(depthIndex)(rowIndices(i) - rowIndex, colIndices(i) - colIndexAtEachRes(depthIndex)))) ;
      }
    }
    rowIndex += nodeToProcess->GetObsInNode().size() ; // The B matrices should have as many rows as observations in the node...
  }
  m_obsOrderForFmat = FmatObsOrder ;

  // Adding in the intercept...
  for (uint i = 0; i < m_dataset.covariateValues.rows(); i++) {
    tripletList.push_back(Triplet(i, 0, 1)) ;
  }
  // Processing other covariate values
  for (uint rowIndex = 0 ; rowIndex < m_dataset.covariateValues.rows() ; rowIndex++) {
    for(uint colIndex = 0 ; colIndex < m_dataset.covariateValues.cols() ; colIndex++) {
      tripletList.push_back(Triplet(rowIndex, colIndex + 1, m_dataset.covariateValues(m_obsOrderForFmat(rowIndex), colIndex))) ;
    }
  }

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
  m_Hmat.resize(m_dataset.covariateValues.rows(), m_numKnots + m_dataset.covariateValues.cols() + 1) ;

  m_Hmat.setFromTriplets(tripletList.begin(), tripletList.end()) ;
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
  int loopIndex = 0 ;

  for (int k = 0; k < m_HmatPred.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_HmatPred, k); it; ++it)
    {
      if (it.col() >= m_fixedEffParameters.size()) { // Covariates never change, only the elements in the F matrix change.
        it.valueRef() = *(m_pointerOffsetForHmatPred.at(loopIndex).matrixLocation->data() + m_pointerOffsetForHmatPred.at(loopIndex).offset) ;
        loopIndex += 1 ;
      }
    }
  }
}

void AugTree::createHmatrixPred() {
  uint numObs = m_predictData.spatialCoords.rows() ;
  uint rowIndex = 0 ;

  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
  int numTips = tipNodes.size() ;

  std::vector<ArrayXi> ancestorIdsVec(numTips) ;

  for (uint i = 0 ; i < tipNodes.size(); i++) {
    ArrayXi idVec = tipNodes.at(i)->GetAncestorIds() ; // Last element is tip node.
    ancestorIdsVec.at(i) = idVec ;
  }

  std::vector<Triplet> tripletList;

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

  ArrayXi colIndexAtEachRes = ArrayXi::Zero(m_M + 1);

  for (auto & node : m_vertexVector) {
    if (node->GetDepth() < m_M) {
      for (uint i = node->GetDepth() + 1; i <= m_M; i++) {
        colIndexAtEachRes(i) += node->GetNumKnots() ;
      }
    }
  } // Convoluted way of getting cumsum(c(0, Num knots at resolutions 0 to M-1))

  colIndexAtEachRes += (m_dataset.covariateValues.cols() + 1) ; // Covariate columns go before, the +1 is for the intercept.
  std::vector<TreeNode *> previousBrickAncestors = tipNodes.at(0)->getAncestors() ;

  for (auto & nodeToProcess : tipNodes) {

    // The idea behind this section of the code is that the column index for producing the section
    // of the H matrix for a tip node at any given depth i should only change when a different ancestor for the
    // tip node is reached. The hierarchical structure of the tree explains this.
    for (uint i = 1 ; i <= m_M; i++) { // The root node is shared by all tips, hence the 1.
      if (previousBrickAncestors.at(i) != nodeToProcess->getAncestors().at(i)) {
        colIndexAtEachRes(i) += previousBrickAncestors.at(i)->GetNumKnots() ;
      }
    }
    previousBrickAncestors = nodeToProcess->getAncestors() ;

    if (nodeToProcess->GetPredIndices().size() == 0 ) {
      continue ;
    }

    for (uint depthIndex = 0; depthIndex <= m_M; depthIndex++) {
      ArrayXi rowIndices = rep(uvec::LinSpaced(nodeToProcess->GetUpred(depthIndex).rows(), 0, nodeToProcess->GetUpred(depthIndex).rows() - 1).array(),
                               nodeToProcess->GetUpred(depthIndex).cols()) + rowIndex  ;
      ArrayXi colIndices = rep_each(uvec::LinSpaced(nodeToProcess->GetUpred(depthIndex).cols(), 0, nodeToProcess->GetUpred(depthIndex).cols() - 1).array(),
                                    nodeToProcess->GetUpred(depthIndex).rows()) + colIndexAtEachRes(depthIndex) ;

      for (uint i = 0; i < rowIndices.size(); i++) {
        tripletList.push_back(Triplet(rowIndices(i), colIndices(i), nodeToProcess->GetUpred(depthIndex)(rowIndices(i) - rowIndex, colIndices(i) - colIndexAtEachRes(depthIndex)))) ;
      }
    }
    rowIndex += nodeToProcess->GetPredIndices().size() ; // The U matrices should have as many rows as prediction locations in the node...
  }

  // Adding in the intercept...
  for (uint i = 0; i < m_predictData.covariateValues.rows(); i++) {
    tripletList.push_back(Triplet(i, 0, 1)) ;
  }
  // Processing other covariate values
  for (uint rowInd = 0 ; rowInd < m_predictData.covariateValues.rows() ; rowInd++) {
    for(uint colInd = 0 ; colInd < m_predictData.covariateValues.cols() ; colInd++) {
      tripletList.push_back(Triplet(rowInd, colInd + 1, m_predictData.covariateValues(m_obsOrderForHpredMat(rowInd), colInd))) ;
    }
  }
  Rcout << "Creating HmatPred! \n" ;
  m_HmatPred.resize(m_predictData.covariateValues.rows(), m_numKnots + m_dataset.covariateValues.cols() + 1) ;
  m_HmatPred.setFromTriplets(tripletList.begin(), tripletList.end()) ;
  Rcout << "Done! \n" ;

  for (auto & tipNode : tipNodes) {
    if (tipNode->GetPredIndices().size() > 0) {
      std::vector<TreeNode *> ancestorsList = tipNode->getAncestors() ;
      uint numObs = tipNode->GetPredIndices().size() ;
      for (uint rowIndex = 0; rowIndex < numObs ; rowIndex++) {
        for (uint depth = 0 ; depth <= m_M ; depth++) {
          uint numKnotsAtDepth = ancestorsList.at(depth)->GetNumKnots() ;
          for (uint colIndex = 0; colIndex < numKnotsAtDepth; colIndex++) {
            pointerOffset elementToAdd = pointerOffset(&(tipNode->GetUmatList().at(depth)), colIndex * numObs + rowIndex) ;
            m_pointerOffsetForHmatPred.push_back(elementToAdd) ;
          }
        }
      }
    }
  }
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

// Am I betraying the concept of unique pointers here?
// Since I create this vector of pairs only to simplify the syntax of the
// final evaluation, maybe it's ok? There might be a better way to go about it though.
void AugTree::ComputeLogPriors() {

  std::vector<std::pair<double, TwoParsProbDist *>> priorCombinations ;

  priorCombinations.push_back(std::make_pair(m_MaternParsSpace.m_rho, m_MaternParsHyperparsRhoSpace.get())) ;
  priorCombinations.push_back(std::make_pair(m_MaternParsSpace.m_smoothness, m_MaternParsHyperparsSmoothnessSpace.get())) ;

  priorCombinations.push_back(std::make_pair(m_MaternParsTime.m_rho, m_MaternParsHyperparsRhoTime.get())) ;
  priorCombinations.push_back(std::make_pair(m_MaternParsTime.m_smoothness, m_MaternParsHyperparsSmoothnessTime.get())) ;

  priorCombinations.push_back(std::make_pair(m_spacetimeScaling, m_MaternParsHyperparsScaling.get())) ;

  priorCombinations.push_back(std::make_pair(m_fixedEffSD, m_ParsHyperparsFixedEffsSD.get())) ;
  priorCombinations.push_back(std::make_pair(m_errorSD, m_ParsHyperparsErrorSD.get())) ;

  double logPrior = 0 ;

  for (auto & i : priorCombinations) {
    double coef = i.first ;
    if (m_normalHyperprior) coef = log(i.first) ; // The normal prior is applicable to log-hyperparameters.
    logPrior += i.second->computeLogDensity(i.first) ;
  }

  m_logPrior = logPrior ;
}

void AugTree::ComputeLogFCandLogCDandDataLL() {

  int n = m_dataset.responseValues.size() ;
  Rcout << "Creating SigmaBetaEtaInv..." << std::endl ;
  if (m_SigmaFEandEtaInv.rows() == 0) {
    CreateSigmaBetaEtaInvMat() ;
  } else {
    UpdateSigmaBetaEtaInvMat() ;
  }
  Rcout << "Done! Creating Hpred..." << std::endl ;
  if (m_processPredictions) {
    ComputeHpred() ;
  }

  // for (auto & i: m_vertexVector) {
  //   if (i->GetDepth() < m_M) { // For internal nodes
  //     i->clearWmatrices() ;
  //   }
  // }
  Rcout << "Done! Computing Chol. SigmaFEeta..." << std::endl ;
  fflush(stdout) ;
  Eigen::SimplicialLDLT<sp_mat> blockChol(m_SigmaFEandEtaInv) ;

  double logDetSigmaKFEinv = blockChol.vectorD().array().log().sum() ;
  Rcout << "Done! Creating H matrix..." << std::endl ;
  if (m_recomputeMRAlogLik) {
    if (m_Hmat.size() == 0) {
      createHmatrix() ;
    } else {
      updateHmatrix() ;
    }
  }
  Rcout << "Done! Scaling responses..." << std::endl ;
  fflush(stdout) ;
  vec responsesReshuffled = elem(m_dataset.responseValues.array(), m_obsOrderForFmat) ;
  mat scaledResponse = std::pow(m_errorSD, -2) * m_Hmat.transpose() * responsesReshuffled ;
  Rcout << "Done! Computing secondTerm..." << std::endl ;
  fflush(stdout)  ;
  sp_mat secondTerm = std::pow(m_errorSD, -2) * (m_Hmat.transpose() * m_Hmat) ;

  Rcout << "Done! analysing sparsity pattern..." << std::endl ;
  fflush(stdout) ;
  if (m_logFullCond == 0) { // This is the first iteration...
    try {
      m_FullCondPrecisionChol.analyzePattern(m_SigmaFEandEtaInv + secondTerm) ;
    } catch(std::bad_alloc& ex) {
      forward_exception_to_r(ex) ;
    }
  }
  Rcout << "Done! Computing Full Cond. Chol..." << std::endl ;
  fflush(stdout) ;
  m_FullCondPrecisionChol.factorize(m_SigmaFEandEtaInv + secondTerm) ; // The sparsity pattern in matToInvert is always the same, notwithstand the hyperparameter values.
  secondTerm.resize(0,0) ; // Making sure it doesn't clog the memory (although I suspect this is not necessary)
  Rcout << "Done! Computing Full Cond. SDs..." << std::endl ;
  fflush(stdout) ;
  ComputeFullCondSDsFE() ;
  Rcout << "Done! Computing updated means..." << std::endl ;
  fflush(stdout) ;
  vec updatedMean = m_FullCondPrecisionChol.solve(scaledResponse) ;
  if(m_FullCondPrecisionChol.info()!=Success) {
    // solving failed
    Rcpp::Rcout << "Solving failed!!!! \n" ;
    throw Rcpp::exception("Leave now... \n") ;
  }

  m_Vstar = updatedMean ; // Assuming there will be an implicit conversion to vec type.

  m_FullCondMean = m_Vstar ;

  vec fixedEffMeans = m_Vstar.head(m_fixedEffParameters.size()) ;
  SetFixedEffParameters(fixedEffMeans) ;
  Rcout << "Done! Computing log-det Q..." << std::endl ;
  fflush(stdout) ;
  double logDetQmat = m_FullCondPrecisionChol.vectorD().array().log().sum() ; // In LDLT, the L matrix has a diagonal with 1s, meaning that its determinant is 1. It follows that the determinant of the original matrix is simply the product of the elements in the D matrix.

  m_logFullCond = 0.5 * logDetQmat ; // Since we arbitrarily evaluate always at the full-conditional mean, the exponential part of the distribution reduces to 0.

  // Computing p(v* | Psi)
  mat logLikTerm ;
  Rcout << "Done! Computing log-lik..." << std::endl ;
  fflush(stdout) ;
  logLikTerm.noalias() = -0.5 * m_Vstar.transpose() * m_SigmaFEandEtaInv.selfadjointView<Upper>() * m_Vstar ;
  Rcout << "Done! Finalising..." << std::endl ;
  fflush(stdout) ;
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

  // m_MaternParsSpace.print("Spatial parameters:") ;
  // m_MaternParsTime.print("Time parameters:") ;
  // Rprintf("Scaling parameter: %.3e \n", m_spacetimeScaling) ;

  ComputeLogPriors() ;

  // if (m_recomputeMRAlogLik) {
    computeWmats() ;
  // }

  ComputeLogFCandLogCDandDataLL() ;

  // Rprintf("Observations log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
  //  m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  m_logJointPsiMarginal = m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond ;
   // Rprintf("Joint value: %.4e \n \n", m_logJointPsiMarginal) ;
}

void AugTree::ComputeHpred() {
  std::vector<TreeNode *> tipNodes = GetTipNodes() ;
  for (auto & i : tipNodes) {
    ArrayXi predictionsInLeaf = i->GetPredIndices() ;
    if (predictionsInLeaf.size() > 0) {
      i->computeUpred(m_MaternParsSpace, m_MaternParsTime, m_spacetimeScaling, m_predictData, m_nuggetSD, m_distMethod) ;
    }
  }

  if (m_HmatPred.rows() == 0) {
    createHmatrixPred() ;
  } else {
    updateHmatrixPred() ;
  }
}

vec AugTree::ComputeEvar() {
  double errorVar = std::pow(m_errorSD, 2) ;
  vec EvarValues = vec::Zero(m_HmatPred.rows()) ;
  int obsIndex = 0 ;

  vec invertedDsqrt = 1/m_FullCondPrecisionChol.vectorD().array().pow(0.5) ;

  #pragma omp parallel for
  for (uint i = 0; i < m_HmatPred.rows(); i++) {
    vec HmatPredRow = m_HmatPred.row(i) ;
    vec solution = HmatPredRow.transpose() * m_FullCondPrecisionChol.solve(HmatPredRow) ;
    EvarValues(i) = solution(0) + errorVar ;
  }
  return EvarValues ;
}

void AugTree::SetMaternPars(const Rcpp::List & MaternPars) {
  List SpaceParas = Rcpp::as<List>(MaternPars["space"]) ;
  List TimeParas = Rcpp::as<List>(MaternPars["time"]) ;
  double scalePara = Rcpp::as<double>(MaternPars["scale"]) ;

  double rhoSpace = Rcpp::as<double>(SpaceParas["rho"]) ;
  double smoothnessSpace = Rcpp::as<double>(SpaceParas["smoothness"]) ;

  double rhoTime = Rcpp::as<double>(TimeParas["rho"]) ;
  double smoothnessTime = Rcpp::as<double>(TimeParas["smoothness"]) ;

  maternVec MaternParsSpace(rhoSpace, smoothnessSpace, 1) ;
  maternVec MaternParsTime(rhoTime, smoothnessTime, 1) ;

  bool test = (fabs(m_spacetimeScaling - scalePara) < epsilon) && (m_MaternParsSpace == MaternParsSpace) && (m_MaternParsTime == MaternParsTime) ;

  if (test) {
    m_recomputeMRAlogLik = false ;
  } else {
    m_recomputeMRAlogLik = true ;
  }
  m_MaternParsSpace = MaternParsSpace ;
  m_MaternParsTime = MaternParsTime ;
  m_spacetimeScaling = scalePara ;
}

void AugTree::ComputeFullCondSDsFE() {
  m_FullCondSDs = vec::Zero(m_Hmat.cols()) ; // We only derive SDs for FEs to save time though.
  #pragma omp parallel for
  for (uint i = 0; i < m_fixedEffParameters.size(); i++) {
    vec solveVector = vec::Zero(m_Hmat.cols()) ;
    solveVector(i) = 1 ;
    vec matrixCol = m_FullCondPrecisionChol.solve(solveVector) ;
    m_FullCondSDs(i) = pow(matrixCol(i), 0.5) ;
  }
}

void AugTree::SetParsHyperpars(const Rcpp::List & MaternParsHyperparsList,
                               const Rcpp::NumericVector & fixedEffsParsHyperpars,
                               const Rcpp::NumericVector & errorParsHyperpars,
                               const bool normalHyperprior) {
  Rcpp::List spaceParas = Rcpp::as<Rcpp::List>(MaternParsHyperparsList["space"]) ;
  Rcpp::List timeParas = Rcpp::as<Rcpp::List>(MaternParsHyperparsList["time"]) ;

  if (normalHyperprior) { // How can this be improved? Very redundant...
    m_MaternParsHyperparsScaling = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(MaternParsHyperparsList["scale"]))) ;
    m_MaternParsHyperparsRhoSpace = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(spaceParas["rho"]))) ;
    m_MaternParsHyperparsSmoothnessSpace = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(spaceParas["smoothness"]))) ;
    m_MaternParsHyperparsRhoTime = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(timeParas["rho"]))) ;
    m_MaternParsHyperparsSmoothnessTime = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(timeParas["smoothness"]))) ;
    m_ParsHyperparsFixedEffsSD = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(fixedEffsParsHyperpars))) ;
    m_ParsHyperparsErrorSD = std::unique_ptr<NormalDist>(new NormalDist(Rcpp::as<vec>(errorParsHyperpars))) ;
  } else {
    m_MaternParsHyperparsScaling = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(MaternParsHyperparsList["scale"]))) ;
    m_MaternParsHyperparsRhoSpace = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(spaceParas["rho"]))) ;
    m_MaternParsHyperparsSmoothnessSpace = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(spaceParas["smoothness"]))) ;
    m_MaternParsHyperparsRhoTime = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(timeParas["rho"]))) ;
    m_MaternParsHyperparsSmoothnessTime = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(timeParas["smoothness"]))) ;
    m_ParsHyperparsFixedEffsSD = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(fixedEffsParsHyperpars))) ;
    m_ParsHyperparsErrorSD = std::unique_ptr<GammaDist>(new GammaDist(Rcpp::as<vec>(errorParsHyperpars))) ;
  }
}

