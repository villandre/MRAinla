#include<cmath>
#include<iostream>
#include<cstdlib>

#include "helper.h"

using namespace arma ;

Spatiotemprange sptimeDistance(const arma::vec & spCoor1, const unsigned int & time1, const arma::vec & spCoor2,
                               const unsigned int & time2) {
  arma::vec diffVec = spCoor1 - spCoor2 ;
  arma::vec scaledVec = arma::pow(diffVec, 2) ;
  double sp = arma::sum(scaledVec) ;
  sp = std::sqrt(sp) ;
  int timeDiff = time2 - time1 ;
  timeDiff = abs(timeDiff) ;
  unsigned int timeDiffUint = (unsigned int) timeDiff ;
  return Spatiotemprange(sp, timeDiffUint) ;
};

arma::sp_mat createSparseMatrix(std::vector<arma::mat *> listOfMatrices) {
  int numMatrices = listOfMatrices.size() ;
  int dimen = 0 ;
  arma::ivec dimvec(numMatrices) ;

  for(unsigned int i = 0; i < numMatrices; i++) {
    dimvec[i] = listOfMatrices.at(i)->n_rows ;
    dimen += dimvec[i] ;
  }

  arma::sp_mat X(dimen, dimen);

  int idx=0;
  for(unsigned int i = 0; i < numMatrices; i++) {
    X(idx, idx, size(*(listOfMatrices.at(i)))) = *(listOfMatrices.at(i)) ;
    idx = idx + dimvec[i] ;
  }

  return X ;
}

std::vector<unsigned int> extractBlockIndices(const arma::sp_mat & symmSparseMatrix) {
  std::vector<unsigned int> blockIndices ;
  blockIndices.push_back(0) ;
  int newPosNextBlock = 0 ;
  int posNextBlock = 0 ;

  while (posNextBlock < symmSparseMatrix.n_cols) {
    int newPosNextBlock = posNextBlock + 1 ; // Ensures that the for loop starts
    printf("Processing block starting at %i ... \n", posNextBlock) ;
    for (int i = posNextBlock; i < newPosNextBlock; i++) {
      arma::vec myCol = arma::vec(symmSparseMatrix.col(i)) ;
      arma::uvec nonZeroElements = arma::find(myCol) ;
      arma::uword lastNonZero = nonZeroElements.tail(1)(0) + 1 ;
      if (lastNonZero > newPosNextBlock) {
        newPosNextBlock = lastNonZero ;
        printf("Moving bound to %i... \n", newPosNextBlock) ;
      }
    }
    blockIndices.push_back(posNextBlock) ;
    posNextBlock = newPosNextBlock ;
  }
  return blockIndices ;
}

// The break is below or to the right of each element. For example, putting rowBreak = 0 and colBreak = 0
// corresponds to creating A11 1 x 1, A12 1 x n_col-1, A21 n_row - 1 x 1, A22 n_row - 1 x n_col - 1
// mat invPDsymmMatrixWithSplit(const sp_mat & pdsymmMatrix, const int rowBreak, const int colBreak) {
//   mat A11 = mat(pdsymmMatrix.submat(0, 0, rowBreak, colBreak)) ;
//   mat A12 = mat(pdsymmMatrix.submat(0, colBreak + 1, rowBreak, pdsymmMatrix.n_cols - 1)) ;
//   mat A21 = mat(pdsymmMatrix.submat(rowBreak + 1, 0, pdsymmMatrix.n_rows, colBreak)) ;
//   mat A22 = mat(pdsymmMatrix.submat(rowBreak + 1, colBreak + 1, pdsymmMatrix.n_rows - 1, pdsymmMatrix.n_cols - 1)) ;
//   // A11 and A22 are square and symmetrical. We check if they are block diagonal.
//   bool blockDiagA11 = blockDiagCheck(A11) ;
//   bool blockDiagA22 = blockDiagCheck(A22) ;
//
//   mat A11inv, A22inv ;
//   if (blockDiagA11) {
//     std::vector<mat> blocksInA11 = extractBlocks(sp_mat(A11)) ;
//     A11inv = invertSymmBlockDiag(blocksInA11) ;
//   } else {
//     A11inv = inv_sympd(A11) ;
//   }
//   if (blockDiagA22) {
//     std::vector<mat> blocksInA22 = extractBlocks(sp_mat(A22)) ;
//     A22inv = invertSymmBlockDiag(blocksInA22) ;
//   } else {
//     A22inv = inv_sympd(A22) ;
//   }
//
//   mat firstInvertedBlock = inv_sympd(A22 - trans(A12) * A11inv * A12) ;
//   mat B11 = A11inv + A11inv * A12 * firstInvertedBlock * trans(A12) * A11inv ;
//   mat B22 = A22inv + A22inv * trans(A12) * inv_sympd(A11 - A12 * A22inv * trans(A12)) * A12 * A22inv ;
//   mat B12t = -A11inv * A12 * firstInvertedBlock ;
//
//   mat Bmatrix = join_rows(join_cols(B11, B12t), join_cols(trans(B12t), B22)) ;
//   return Bmatrix ;
// }

mat invertSymmBlockDiag(const sp_mat & blockMatrix, const uvec & blockIndices) {
  int numRows = blockMatrix.n_rows ;
  mat inverted(numRows, numRows, fill::zeros) ;

  for (unsigned int i = 0; i < (blockIndices.size()-1); i++) {
    unsigned int blockSize = blockIndices.at(i+1) - blockIndices.at(i) ;
    inverted.submat(i, i, i+blockSize-1, i+blockSize-1) =
      inv_sympd(mat(blockMatrix.submat(i, i, i+blockSize-1, i+blockSize-1))) ;
  }
  return inverted ;
}

// bool blockDiagCheck(const mat & matrixToCheck) {
//   int newPosNextBlock = 1 ;
//
//   for (int i = 0; i < newPosNextBlock; i++) {
//     arma::vec myCol = arma::vec(matrixToCheck.col(i)) ;
//     arma::uvec nonZeroElements = arma::find(myCol) ;
//     arma::uword lastNonZero = nonZeroElements.tail(1)(0) + 1 ;
//     if (lastNonZero > newPosNextBlock) {
//       newPosNextBlock = lastNonZero ;
//     }
//   }
//   bool testResult = newPosNextBlock >= matrixToCheck.n_rows ;
//   return testResult ;
// }



