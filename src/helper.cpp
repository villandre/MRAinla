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
typedef Eigen::VectorXd vec ;
typedef Eigen::MatrixXd mat ;
typedef Eigen::SparseMatrix<double> sp_mat ;
typedef Eigen::VectorXi uvec ;
typedef Eigen::MatrixXi umat ;
typedef Eigen::Triplet<double> Triplet;

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

  std::vector<Triplet> tripletList;

  int offset = 0 ;
  for (auto & aMatrix : listOfMatrices) {
    for (uint j = 0; j < aMatrix->cols(); j++) {
      for (uint i = 0; i < aMatrix->rows() ; i++) {
        tripletList.push_back(Triplet(i,j,(*aMatrix(i + offset,j + offset)))) ;
      }
    }
    offset += aMatrix->rows() ;
  }
  sp_mat X(offset, offset);
  X.setFromTriplets(tripletList.begin(), tripletList.end()) ;

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
      bool nonZeroCheck = false ;
      int lastNonZero = myCol.size()-1 ;
      for (auto & i : myCol.reverse()) {
        if (i > 0)
      }

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
    logValue += (-log(sd(i)) -
      0.5 * pow((x(i) - mu(i)), 2)/pow(sd(i), 2)) ;
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
    unsigned int matSize = blockIndices(i+1) - blockIndices(i) ;
    log_det(value, sign, mat(blockMatrix(blockIndices(i), blockIndices(i), size(matSize, matSize)))) ;
    if (sign < 0) {
      throw Rcpp::exception("Error logDetBlockMatrix: Log determinant sign should be positive. \n") ;
    }
    logDeterminant += value ;
  }
  return logDeterminant ;
}

Eigen::Matrix<bool, Eigen::Dynamic, 1> find(const Eigen::Ref<const Eigen::Matrix<bool, Eigen::Dynamic, 1>> & logicalVector) {
  Eigen::VectorXi outputVec(logicalVector.size()) ;
  uint index = 0 ;
  for (uint i = 0 ; i < logicalVector.size(); i++) {
    if (logicalVector(i)) {
      outputVec(index) = i ;
      index += 1 ;
    }
  }
  return outputVec.segment(0, index) ;
}

Eigen::Matrix<bool, Eigen::Dynamic, 1> operator==(const Eigen::Ref<const Eigen::VectorXi> & EigenVec, const uint constant) {
  Eigen::Matrix<bool, Eigen::Dynamic, 1> container(EigenVec.size()) ;
  for (uint i = 0; i < EigenVec.size(); i++) {
    container(i) = EigenVec(i) == constant ;
  }
  return container ;
}
