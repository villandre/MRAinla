// [[Rcpp::depends(BH)]]

#include<math.h>
#include<iostream>
#include<stdlib.h>
#include<gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/bessel.hpp>

#include "helper.h"

using namespace arma ;
using namespace boost ;
using namespace std ;
using namespace math ;

Spatiotemprange sptimeDistance(const arma::vec & spCoor1, const double & time1, const arma::vec & spCoor2,
                               const double & time2) {
  arma::vec diffVec = spCoor1 - spCoor2 ;
  arma::vec scaledVec = arma::pow(diffVec, 2) ;
  double sp = arma::sum(scaledVec) ;
  sp = std::sqrt(sp) ;
  double timeDiff = abs(time2 - time1) ;
  // printf("Time coordinates and time difference: %.4e  %.4e %.4e \n", time1, time2, timeDiff) ;
  return Spatiotemprange(sp, timeDiff) ;
};

arma::sp_mat createBlockMatrix(std::vector<arma::mat *> listOfMatrices) {
  int numMatrices = listOfMatrices.size() ;
  int dimen = 0 ;
  for (auto & i : listOfMatrices) dimen += i->n_rows ;
  arma::ivec dimvec(numMatrices) ;
  arma::umat posMatrix(2, 0) ;
  int rowIndex = 0 ;
  int colIndex = 0 ;
  vec concatenatedValues ;
  for (auto & matrix : listOfMatrices) {
    uvec rowIndices = rep(regspace<uvec>(0, matrix->n_rows - 1), matrix->n_cols) + rowIndex  ;
    uvec colIndices = rep_each(regspace<uvec>(0, matrix->n_cols - 1), matrix->n_rows) + colIndex ;
    posMatrix = join_rows(posMatrix, trans(join_rows(rowIndices, colIndices))) ;
    concatenatedValues = join_cols(concatenatedValues, arma::vectorise(*matrix)) ;
    rowIndex += matrix->n_rows ;
    colIndex += matrix->n_cols ;
  }

  arma::sp_mat X(posMatrix, concatenatedValues, false) ;

  // arma::sp_mat X(dimen, dimen);

  // int idx=0;
  // for(int i = 0; i < numMatrices; i++) {
  //   sp_mat dereferencedMatrix = conv_to<sp_mat>::from(*(listOfMatrices.at(i))) ;
  //   X(idx, idx, arma::size(dereferencedMatrix.n_rows, dereferencedMatrix.n_cols)) = dereferencedMatrix ;
  //   idx += dereferencedMatrix.n_rows ;
  // }

  return X ;
}

arma::uvec extractBlockIndices(const arma::sp_mat & symmSparseMatrix) {
  std::vector<unsigned int> blockIndices ;
  blockIndices.push_back(0) ;

  int posNextBlock = 0 ;
  int newPosNextBlock = 0 ;

  while (posNextBlock < symmSparseMatrix.n_cols) {
    // printf("Processing block starting at %i ... \n", posNextBlock) ;
    newPosNextBlock = posNextBlock + 1 ;
    for (int i = posNextBlock; i < newPosNextBlock; i++) {
      arma::vec myCol = arma::vec(symmSparseMatrix.col(i)) ;
      arma::uvec nonZeroElements = arma::find(myCol) ;
      arma::uword lastNonZero = nonZeroElements.tail(1)(0) ;
      if ((lastNonZero + 1) > newPosNextBlock) {
        newPosNextBlock = lastNonZero + 1;
        // printf("Moving bound to %i... \n", newPosNextBlock) ;
      }
    }
    posNextBlock = newPosNextBlock ;
    blockIndices.push_back(posNextBlock) ;
  }
  return conv_to<arma::uvec>::from(blockIndices) ;
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

sp_mat invertSymmBlockDiag(const sp_mat & blockMatrix, const uvec & blockIndices) {
  int numRows = blockMatrix.n_rows ;

  unsigned int diagElement = 0 ;
  unsigned int blockSize ;
  arma::umat posMatrix(2, 0) ;
  vec concatenatedValues ;
  uvec blockSizes = diff(blockIndices) ;
  int index = 0 ;

  for (uint i = 0 ; i < blockSizes.size() ; i++) {
    umat posMat = join_rows(rep(regspace<uvec>(0, blockSizes.at(i) - 1), blockSizes.at(i)) + index,
              rep_each(regspace<uvec>(0, blockSizes.at(i) - 1), blockSizes.at(i)) + index) ;
    posMatrix = join_rows(posMatrix, trans(posMat)) ;
    index += blockSizes.at(i) ;
  }

  // for (unsigned int i = 0; i < (blockIndices.size()-1); i++) {
  //   blockSize = blockIndices.at(i+1) - blockIndices.at(i) ;
  //   inverted(diagElement, diagElement, size(blockSize, blockSize)) =
  //     inv_sympd(mat(blockMatrix(diagElement, diagElement, size(blockSize, blockSize)))) ;
  //   diagElement += blockSize ;
  // }

  for (unsigned int i = 0; i < (blockIndices.size()-1); i++) {
    blockSize = blockIndices.at(i+1) - blockIndices.at(i) ;
    concatenatedValues = join_cols(concatenatedValues,
                                   vectorise(inv_sympd(mat(blockMatrix(diagElement, diagElement, size(blockSize, blockSize)))))) ;
    diagElement += blockSize ;
  }
  sp_mat inverted(posMatrix, concatenatedValues, false) ;

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

double logNormPDF(const arma::vec & x, const arma::vec & mu, const arma::vec & sd) {
  double logValue = 0;
  for (unsigned int i = 0 ; i < x.size() ; i++) {
    logValue += (-log(sd.at(i)) -
      0.5 * pow((x.at(i) - mu.at(i)), 2)/pow(sd.at(i), 2)) ;
  }
  logValue += (-0.5*x.size()*(log(2) + log(PI))) ;
  return logValue ;
}

// The notation is borrowed from the Wikipedia entry.

double maternCov(const double & distance, const double & rho,
                    const double & smoothness, const double & scale, const double & nugget) {
  double maternValue ;

  if (distance == 0) {
    maternValue = pow(scale, 2) + nugget ;
  } else {
    double base = pow(2 * smoothness, 0.5) * distance / rho ;
    double bessel = boost::math::cyl_bessel_k(smoothness, base) ;
    maternValue = pow(scale, 2) * pow(2, 1 - smoothness) / gsl_sf_gamma(smoothness) *
      pow(base, smoothness) * bessel ;
  }

  return maternValue ;
}

double logDetBlockMatrix(const arma::sp_mat & blockMatrix, const arma::uvec & blockIndices) {
  // uvec blockIndices = extractBlockIndices(blockMatrix) ;
  double logDeterminant = 0 ;
  for (unsigned int i = 0 ; i < (blockIndices.size() - 1) ; i++) {
    double value = 0 ;
    double sign = 0 ;
    unsigned int matSize = blockIndices.at(i+1) - blockIndices.at(i) ;
    log_det(value, sign, mat(blockMatrix(blockIndices.at(i), blockIndices.at(i), size(matSize, matSize)))) ;
    if (sign < 0) {
      throw Rcpp::exception("Error logDetBlockMatrix: Log determinant sign should be positive. \n") ;
    }
    logDeterminant += value ;
  }
  return logDeterminant ;
}
