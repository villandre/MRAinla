// [[Rcpp::depends(BH)]]

#include<math.h>
#include<iostream>
#include<stdlib.h>
#include<gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/bessel.hpp>

#include "helper.h"

using namespace boost ;
using namespace std ;
using namespace math ;
typedef Eigen::VectorXf vec ;
typedef Eigen::MatrixXf mat ;
typedef Eigen::SparseMatrix<float> sp_mat ;
typedef Eigen::VectorXi uvec ;
typedef Eigen::MatrixXi umat ;

Spatiotemprange sptimeDistance(const vec & spCoor1, const double & time1, const vec & spCoor2,
                               const double & time2) {
  vec diffVec = spCoor1 - spCoor2 ;
  vec scaledVec = diffVec.pow(2) ;
  double sp = scaledVec.sum() ;
  sp = std::sqrt(sp) ;
  double timeDiff = abs(time2 - time1) ;
  // printf("Time coordinates and time difference: %.4e  %.4e %.4e \n", time1, time2, timeDiff) ;
  return Spatiotemprange(sp, timeDiff) ;
};

// Pretty slow. Should not be called too often.

sp_mat createBlockMatrix(std::vector<mat *> listOfMatrices) {
  uint numRows = 0 ;
  uint numCols = 0 ;
  for (auto & i : listOfMatrices) {
    numRows += i->n_rows ;
    numCols += i->n_cols ;
  }
  sp_mat X(numRows, numCols);

  int idxRows = 0;
  int idxCols = 0 ;

  for(int i = 0; i < listOfMatrices.size(); i++) {
    sp_mat dereferencedMatrix = conv_to<sp_mat>::from(*(listOfMatrices.at(i))) ;
    X(idxRows, idxCols, arma::size(dereferencedMatrix.n_rows, dereferencedMatrix.n_cols)) = dereferencedMatrix ;
    idxRows += dereferencedMatrix.n_rows ;
    idxCols += dereferencedMatrix.n_cols ;
  }

  return X ;
}

uvec extractBlockIndices(const sp_mat & symmSparseMatrix) {
  std::vector<unsigned int> blockIndices ;
  blockIndices.push_back(0) ;

  int posNextBlock = 0 ;
  int newPosNextBlock = 0 ;

  while (posNextBlock < symmSparseMatrix.n_cols) {
    // printf("Processing block starting at %i ... \n", posNextBlock) ;
    newPosNextBlock = posNextBlock + 1 ;
    for (int i = posNextBlock; i < newPosNextBlock; i++) {
      vec myCol = vec(symmSparseMatrix.col(i)) ;
      uvec nonZeroElements = arma::find(myCol) ;
      int lastNonZero = nonZeroElements.tail(1)(0) ;
      if ((lastNonZero + 1) > newPosNextBlock) {
        newPosNextBlock = lastNonZero + 1;
        // printf("Moving bound to %i... \n", newPosNextBlock) ;
      }
    }
    posNextBlock = newPosNextBlock ;
    blockIndices.push_back(posNextBlock) ;
  }
  return conv_to<uvec>::from(blockIndices) ;
}

double logNormPDF(const vec & x, const vec & mu, const vec & sd) {
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
    if (smoothness >= 1e6) {
      maternValue = pow(scale, 2) * exp(-pow(distance, 2)/(2 * pow(rho, 2))) ;
    } else if (smoothness == 0.5) {
      maternValue = pow(scale, 2) * exp(-distance / rho) ;
    } else if (smoothness == 1.5) {
      maternValue = pow(scale, 2) *
        (1 + sqrt(3) * distance / rho) *
        exp(-sqrt(3) * distance / rho) ;
    } else if (smoothness == 2.5) {
      maternValue = pow(scale, 2) *
        (1 + sqrt(5) * distance/rho + 5 * pow(distance, 2)/(3 * pow(rho, 2))) *
        exp(-sqrt(5) * distance/rho) ;
    } else {
      double base = pow(2 * smoothness, 0.5) * distance / rho ;
      double bessel = boost::math::cyl_bessel_k(smoothness, base) ;
      maternValue = pow(scale, 2) * pow(2, 1 - smoothness) / gsl_sf_gamma(smoothness) *
        pow(base, smoothness) * bessel ;
    }
  }
  return maternValue ;
}

double logDetBlockMatrix(const sp_mat & blockMatrix, const uvec & blockIndices) {
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
