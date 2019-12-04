#include "Forest.h"

using namespace Rcpp ;
using namespace MRAinla ;
using namespace Eigen ;
using namespace std ;

Forest::Forest(uint & Mlon,
                 uint & Mlat,
                 Array2d & lonRange,
                 Array2d & latRange,
                 vec & responses,
                 ArrayXXd & obsSp,
                 ArrayXd & obsTime,
                 ArrayXXd & predCovariates,
                 ArrayXXd & predSp,
                 ArrayXd & predTime,
                 unsigned long int & seed,
                 ArrayXXd & covariates,
                 const unsigned int numKnotsRes0,
                 double J,
                 const string & distMethod)
  : m_distMethod(distMethod)
{
  m_GammaParasSet = false ;
  m_dataset = inputdata(responses, obsSp, obsTime, covariates) ;
  m_mapDimensions = spaceDimensions(lonRange, latRange) ;
  std::random_device rd ;
  std::mt19937_64 generator(rd()) ;
  generator.seed(seed) ;
  m_randomNumGenerator = generator ;
  SetPredictData(predSp, predTime, predCovariates) ;
  m_assignedPredToKnot = Array<bool, Dynamic, 1>(m_predictData.timeCoords.size()) ;
  m_assignedPredToKnot.segment(0, m_assignedPredToKnot.size()) = false ;

  m_fixedEffParameters = Eigen::VectorXd::Zero(m_dataset.covariateValues.cols() + 1);

  Rcout << "Building the grid..." << std::endl ;
  std::vector<double> timeAsVector;
  timeAsVector.resize(obsTime.size());
  ArrayXd::Map(&timeAsVector[0], obsTime.size()) = obsTime ;
  std::sort(timeAsVector.begin(), timeAsVector.end()) ;
  timeAsVector.erase(unique( timeAsVector.begin(), timeAsVector.end() ),
                     timeAsVector.end() );
  std::vector<AugTree *> MRAgridVector ;

  for (auto & i : timeAsVector) {
    m_treeVector.push_back(AugTree(Mlon, Mlat, m_mapDimensions, i, m_dataset, m_predictData, numKnotsRes0, J, m_distMethod, m_randomNumGenerator, m_assignedPredToKnot)) ;
  }

  for (auto & tree : m_treeVector) {
    m_numKnots += tree.GetNumKnots() ;
  }

  m_FullCondMean = Eigen::VectorXd::Zero(m_numKnots + m_fixedEffParameters.size() + m_uniqueTimeValues.size()) ;
  m_FullCondSDs = Eigen::VectorXd::Zero(m_numKnots + m_fixedEffParameters.size() + m_uniqueTimeValues.size()) ;
  m_Vstar = Eigen::VectorXd::Zero(m_numKnots +  m_fixedEffParameters.size() + m_uniqueTimeValues.size()) ;
  m_Vstar.segment(0, m_FEmu.size()) = m_FEmu ;
}

void Forest::CreateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList ;
  for (auto & tree : m_treeVector) {
    std::vector<mat *> KinvVec = tree.getKmatricesInversePointers() ;
    FEinvAndKinvMatrixList.insert( FEinvAndKinvMatrixList.end(), KinvVec.begin(), KinvVec.end() );
  }
  mat timeMat = createTimePrecisionMatrix(m_timeCovPara, m_uniqueTimeValues) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &timeMat) ;
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

void Forest::UpdateSigmaBetaEtaInvMat() {
  std::vector<mat *> FEinvAndKinvMatrixList ;
  for (auto & tree : m_treeVector) {
    std::vector<mat *> KinvVec = tree.getKmatricesInversePointers() ;
    FEinvAndKinvMatrixList.insert( FEinvAndKinvMatrixList.end(), KinvVec.begin(), KinvVec.end() );
  }
  mat timeMat = createTimePrecisionMatrix(m_timeCovPara, m_uniqueTimeValues) ;
  FEinvAndKinvMatrixList.insert(FEinvAndKinvMatrixList.begin(), &timeMat) ;
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

void Forest::createHmatrix() {
  // int numObs = m_dataset.spatialCoords.rows() ;
  uint bigColOffset = 0 ;
  ArrayXi FmatObsOrder = ArrayXi::Zero(m_dataset.timeCoords.size()) ;
  int rowIndex = 0 ;
  // In this new scheme, we have columns representing time ;
  m_Hmat.resize(m_dataset.covariateValues.rows(), m_numKnots + m_uniqueTimeValues.size() + m_dataset.covariateValues.cols() + 1) ;
  std::vector<Triplet> tripletList ;

  for (auto & tree : m_treeVector) {
    std::vector<TreeNode *> tipNodes = tree.GetTipNodes() ;
    int numObs = tree.GetNumObs() ;
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

    ArrayXi colIndexAtEachRes = ArrayXi::Zero(tree.GetM() + 1);

    for (auto & node : tree.GetVertexVector()) {
      if (node->GetDepth() < tree.GetM()) {
        for (uint i = node->GetDepth() + 1; i <= tree.GetM(); i++) {
          colIndexAtEachRes(i) += node->GetNumKnots() ;
        }
      }
    } // Convoluted way of getting cumsum(c(0, Num knots at resolutions 0 to M-1))

    colIndexAtEachRes += (m_dataset.covariateValues.cols() + 1) + bigColOffset + m_uniqueTimeValues.size() ; // Covariate columns go before, the +1 is for the intercept.
    std::vector<TreeNode *> previousBrickAncestors = tipNodes.at(0)->getAncestors() ;

    for (auto & nodeToProcess : tipNodes) {
      if (nodeToProcess->GetObsInNode().size() == 0 ) {
        continue ;
      }

      // The idea behind this section of the code is that the column index for producing the section
      // of the H matrix for a tip node at any given depth i should only change when a different ancestor for the
      // tip node is reached. The hierarchical structure of the tree explains this.
      for (uint i = 1 ; i <= tree.GetM(); i++) { // The root node is shared by all tips, hence the 1.
        if (previousBrickAncestors.at(i) != nodeToProcess->getAncestors().at(i)) {
          colIndexAtEachRes(i) += previousBrickAncestors.at(i)->GetNumKnots() ;
        }
      }
      previousBrickAncestors = nodeToProcess->getAncestors() ;

      FmatObsOrder.segment(rowIndex, nodeToProcess->GetObsInNode().size()) = nodeToProcess->GetObsInNode() ;

      for (uint depthIndex = 0; depthIndex <= tree.GetM(); depthIndex++) {
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

    // The idea behind this is to associate a memory location in the W matrices and a cell
    // in the H matrix. The memory location we get from calling .data() on a given W
    // matrix changes when the W values are updated, even when W is allocated on the
    // heap, i.e. with 'new', but the memory location of the vector element is constant
    // once the vector is fully populated. To get access to the value, all we need to do
    // is invoke *(pointer->data() + offset).
    // The odd nesting of the loops is due to the fact that m_Hmat is row-major.
    // I made m_Hmat row-major to make the updating process easier, as all values on any
    // given row are associated with the same tip.
    // Don't mind the four nested loops: all it does is traverse all the elements once.

    for (auto & tipNode : tipNodes) {
      uint numObs = tipNode->GetObsInNode().size() ;
      for (uint rowIndex = 0; rowIndex < numObs ; rowIndex++) {
        for (uint depth = 0 ; depth <= tree.GetM() ; depth++) {
          uint numKnotsAtDepth = tipNode->GetWlist().at(depth).cols() ;
          for (uint colIndex = 0; colIndex < numKnotsAtDepth; colIndex++) {
            pointerOffset elementToAdd = pointerOffset(&(tipNode->GetWlist().at(depth)), colIndex * numObs + rowIndex) ;
            m_pointerOffsetForHmat.push_back(elementToAdd) ;
          }
        }
      }
    }
    bigColOffset += tree.GetNumKnots() ;
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
  // Processing time
  for (uint i = 0 ; i < m_dataset.responseValues.size(); i++) {
    double timeValue = m_dataset.timeCoords(m_obsOrderForFmat(i)) ;
    Array<bool, Dynamic, 1> forFindFun = (timeValue == m_uniqueTimeValues) ;
    ArrayXi colIndex = find(forFindFun) ;
    tripletList.push_back(Triplet(i, colIndex(0) + m_dataset.covariateValues.cols() + 1, timeValue)) ;
  }
  m_Hmat.setFromTriplets(tripletList.begin(), tripletList.end()) ;
}

void Forest::updateHmatrix() {
  int loopIndex = 0 ;

  for (int k=0; k<m_Hmat.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_Hmat, k); it; ++it)
    {
      if (it.col() >= m_fixedEffParameters.size() + m_uniqueTimeValues.size()) { // Covariates and time values never change, only the elements in the F matrices change.
        it.valueRef() = *(m_pointerOffsetForHmat.at(loopIndex).matrixLocation->data() + m_pointerOffsetForHmat.at(loopIndex).offset) ;
        loopIndex += 1 ;
      }
    }
  }
}

void Forest::updateHmatrixPred() {
  int loopIndex = 0 ;

  for (int k = 0; k < m_HmatPred.outerSize(); ++k) {
    for (sp_mat::InnerIterator it(m_HmatPred, k); it; ++it)
    {
      if (it.col() >= m_fixedEffParameters.size() + m_uniqueTimeValues.size()) { // Covariates and time values never change, only the elements in the F matrices change.
        it.valueRef() = *(m_pointerOffsetForHmatPred.at(loopIndex).matrixLocation->data() + m_pointerOffsetForHmatPred.at(loopIndex).offset) ;
        loopIndex += 1 ;
      }
    }
  }
}

void Forest::createHmatrixPred() {
  uint numObs = m_predictData.spatialCoords.rows() ;
  uint rowIndex = 0 ;
  uint bigColOffset = 0 ;
  std::vector<Triplet> tripletList;
  m_obsOrderForFpredMat = ArrayXi::Zero(m_predictData.timeCoords.size()) ;

  for (auto & tree : m_treeVector) {
    std::vector<TreeNode *> tipNodes = tree.GetTipNodes() ;
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

    uint elementIndex = 0 ;
    for (auto & i : tipNodes) {
      uint numPredObsInNode = i->GetPredIndices().size() ;
      if (numPredObsInNode > 0) {
        m_obsOrderForFpredMat.segment(elementIndex, i->GetPredIndices().size()) = i->GetPredIndices() ;
        elementIndex += i->GetPredIndices().size() ;
      }
    }

    ArrayXi colIndexAtEachRes = ArrayXi::Zero(tree.GetM() + 1);

    for (auto & node : tree.GetVertexVector()) {
      if (node->GetDepth() < tree.GetM()) {
        for (uint i = node->GetDepth() + 1; i <= tree.GetM(); i++) {
          colIndexAtEachRes(i) += node->GetNumKnots() ;
        }
      }
    } // Convoluted way of getting cumsum(c(0, Num knots at resolutions 0 to M-1))

    colIndexAtEachRes += (m_dataset.covariateValues.cols() + 1 + m_uniqueTimeValues.size() + bigColOffset) ; // Covariate and time columns go before, the +1 is for the intercept.
    std::vector<TreeNode *> previousBrickAncestors = tipNodes.at(0)->getAncestors() ;

    for (auto & nodeToProcess : tipNodes) {

      if (nodeToProcess->GetPredIndices().size() == 0 ) {
        continue ;
      }

      // The idea behind this section of the code is that the column index for producing the section
      // of the H matrix for a tip node at any given depth i should only change when a different ancestor for the
      // tip node is reached. The hierarchical structure of the tree explains this.
      for (uint i = 1 ; i <= tree.GetM(); i++) { // The root node is shared by all tips, hence the 1.
        if (previousBrickAncestors.at(i) != nodeToProcess->getAncestors().at(i)) {
          colIndexAtEachRes(i) += previousBrickAncestors.at(i)->GetNumKnots() ;
        }
      }
      previousBrickAncestors = nodeToProcess->getAncestors() ;

      for (uint depthIndex = 0; depthIndex <= tree.GetM(); depthIndex++) {
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

    for (auto & tipNode : tipNodes) {
      if (tipNode->GetPredIndices().size() > 0) {
        std::vector<TreeNode *> ancestorsList = tipNode->getAncestors() ;
        uint numObs = tipNode->GetPredIndices().size() ;
        for (uint rowIndex = 0; rowIndex < numObs ; rowIndex++) {
          for (uint depth = 0 ; depth <= tree.GetM() ; depth++) {
            uint numKnotsAtDepth = ancestorsList.at(depth)->GetNumKnots() ;
            for (uint colIndex = 0; colIndex < numKnotsAtDepth; colIndex++) {
              pointerOffset elementToAdd = pointerOffset(&(tipNode->GetUmatList().at(depth)), colIndex * numObs + rowIndex) ;
              m_pointerOffsetForHmatPred.push_back(elementToAdd) ;
            }
          }
        }
      }
    }
    bigColOffset += tree.GetNumKnots() ;
  }
  // Adding in the intercept...
  for (uint i = 0; i < m_predictData.covariateValues.rows(); i++) {
    tripletList.push_back(Triplet(i, 0, 1)) ;
  }
  // Processing other covariate values
  for (uint rowInd = 0 ; rowInd < m_predictData.covariateValues.rows() ; rowInd++) {
    for(uint colInd = 0 ; colInd < m_predictData.covariateValues.cols() ; colInd++) {
      tripletList.push_back(Triplet(rowInd, colInd + 1, m_predictData.covariateValues(m_obsOrderForFpredMat(rowInd), colInd))) ;
    }
  }

  // Processing time
  for (uint i = 0 ; i < m_predictData.responseValues.size(); i++) {
    double timeValue = m_predictData.timeCoords(m_obsOrderForFpredMat(i)) ;
    Array<bool, Dynamic, 1> forFindFun = (timeValue == m_uniqueTimeValues) ;
    ArrayXi colIndex = find(forFindFun) ;
    tripletList.push_back(Triplet(i, colIndex(0) + m_dataset.covariateValues.cols() + 1, timeValue)) ;
  }

  m_HmatPred.resize(m_predictData.covariateValues.rows(), m_numKnots + m_predictData.covariateValues.cols() + 1) ;
  m_HmatPred.setFromTriplets(tripletList.begin(), tripletList.end()) ;
}

// For now, we assume that all hyperpriors have an inverse gamma distribution with the same parameters.

void Forest::ComputeLogPriors() {

  std::vector<std::pair<double, GammaHyperParas>> priorCombinations ;

  priorCombinations.push_back(std::make_pair(m_MRAcovParasSpace.m_rho, m_maternParasGammaAlphaBetaSpace.m_rho)) ;
  priorCombinations.push_back(std::make_pair(m_MRAcovParasSpace.m_smoothness, m_maternParasGammaAlphaBetaSpace.m_smoothness)) ;
  priorCombinations.push_back(std::make_pair(m_MRAcovParasSpace.m_scale, m_maternParasGammaAlphaBetaSpace.m_scale)) ;

  priorCombinations.push_back(std::make_pair(m_fixedEffSD, m_fixedEffGammaAlphaBeta)) ;
  priorCombinations.push_back(std::make_pair(m_errorSD, m_errorGammaAlphaBeta)) ;

  priorCombinations.push_back(std::make_pair(m_timeCovPara, m_timeCovParasGammaAlphaBeta)) ;

  double logPrior = 0 ;

  for (auto & i : priorCombinations) {
    logPrior += (i.second.m_alpha - 1) * log(i.first) - i.second.m_beta * i.first ;
  }

  m_logPrior = logPrior ;
}

void Forest::ComputeLogFCandLogCDandDataLL() {

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
  mat scaledResponse = std::pow(m_errorSD, -2) * responsesReshuffled.transpose() * m_Hmat ;
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
  vec updatedMean = m_FullCondPrecisionChol.solve(scaledResponse.transpose()) ;
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

void Forest::ComputeLogJointPsiMarginal() {

  // m_MRAcovParasSpace.print("Spatial parameters:") ;
  // m_MRAcovParasTime.print("Time parameters:") ;

  ComputeLogPriors() ;

  // if (m_recomputeMRAlogLik) {
  for (auto & tree : m_treeVector) {
    tree.computeWmats(m_MRAcovParasSpace, m_spaceNuggetSD, m_distMethod) ;
  }
  // }

  ComputeLogFCandLogCDandDataLL() ;

  // Rprintf("Observations log-lik: %.4e \n Log-prior: %.4e \n Log-Cond. dist.: %.4e \n Log-full cond.: %.4e \n \n \n",
  //  m_globalLogLik, m_logPrior, m_logCondDist, m_logFullCond) ;
  m_logJointPsiMarginal = m_globalLogLik + m_logPrior + m_logCondDist - m_logFullCond ;
  // Rprintf("Joint value: %.4e \n \n", m_logJointPsiMarginal) ;
}

void Forest::ComputeHpred() {
  for (auto & tree : m_treeVector) {
    std::vector<TreeNode *> tipNodes = tree.GetTipNodes() ;
    for (auto & i : tipNodes) {
      ArrayXi predictionsInLeaf = i->GetPredIndices() ;
      if (predictionsInLeaf.size() > 0) {
        i->computeUpred(m_MRAcovParasSpace, m_predictData, m_spaceNuggetSD, m_distMethod) ;
      }
    }
  }

  if (m_HmatPred.rows() == 0) {
    createHmatrixPred() ;
  } else {
    updateHmatrixPred() ;
  }
}

vec Forest::ComputeEvar() {
  double errorVar = std::pow(m_errorSD, 2) ;
  vec EvarValues = vec::Zero(m_HmatPred.rows()) ;
  int obsIndex = 0 ;

  vec invertedDsqrt = 1/m_FullCondPrecisionChol.vectorD().array().pow(0.5) ;

#pragma omp parallel for
  for (uint i = 0; i < m_HmatPred.rows(); i++) {
    vec HmatPredRow = m_HmatPred.row(i) ;
    vec solution = HmatPredRow.transpose() * m_FullCondPrecisionChol.solve(HmatPredRow) ;
    EvarValues(i) = solution(0) ;
  }
  return EvarValues ;
}

void Forest::SetMRAcovParas(const Rcpp::List & MRAcovParas) {
  List SpaceParas = Rcpp::as<List>(MRAcovParas["space"]) ;

  double rhoSpace = Rcpp::as<double>(SpaceParas["rho"]) ;
  double smoothnessSpace = Rcpp::as<double>(SpaceParas["smoothness"]) ;
  double scaleSpace = Rcpp::as<double>(SpaceParas["scale"]) ;

  maternVec MRAcovParasSpace(rhoSpace, smoothnessSpace, scaleSpace) ;

  bool test = m_MRAcovParasSpace == MRAcovParasSpace ;

  if (test) {
    m_recomputeMRAlogLik = false ;
  } else {
    m_recomputeMRAlogLik = true ;
  }
  m_MRAcovParasSpace = MRAcovParasSpace ;
}

void Forest::SetMRAcovParasGammaAlphaBeta(const Rcpp::List & MRAcovParasList) {
  Rcpp::List spaceParas = Rcpp::as<List>(MRAcovParasList["space"]) ;
  m_maternParasGammaAlphaBetaSpace = maternGammaPriorParas(GammaHyperParas(Rcpp::as<vec>(spaceParas["rho"])),
                                                            GammaHyperParas(Rcpp::as<vec>(spaceParas["smoothness"])),
                                                            GammaHyperParas(Rcpp::as<vec>(spaceParas["scale"]))) ;
}

void Forest::SetTimeCovParaGammaAlphaBeta(const Rcpp::List & timeCovParasList) {
  Rcpp::List timeParas = Rcpp::as<List>(timeCovParasList) ;
  m_timeCovParasGammaAlphaBeta = GammaHyperParas(Rcpp::as<vec>(timeParas["scale"])) ; // Only a scale parameter for now.
}

void Forest::ComputeFullCondSDsFE() {
  mat identityForSolve = mat::Identity(m_FullCondPrecisionChol.vectorD().size(), m_fixedEffParameters.size()) ;
  m_FullCondPrecisionChol.matrixL().solveInPlace(identityForSolve) ;

  vec invertedDiag = m_FullCondPrecisionChol.vectorD().array().pow(-1).matrix() ;

  m_FullCondSDs = vec::Zero(m_FullCondPrecisionChol.vectorD().size()) ;
  identityForSolve.noalias() = identityForSolve.array().pow(2).matrix() ;
  vec varValues = identityForSolve * invertedDiag ;
  m_FullCondSDs.segment(0, m_fixedEffParameters.size()) = varValues.array().pow(0.5) ;
}
