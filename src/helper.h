#include<RcppEigen.h>
#include<string.h>
#include <random>

#ifndef HELPER_H
#define HELPER_H

struct Spatiotemprange{
  double sp ;
  double time ;

  Spatiotemprange(double & sp, double & time) : sp(sp), time(time) { } ;
  Spatiotemprange() { } ;
  void print() const {
    Rprintf("Spatial, temporal range: %.3e, %.3e \n", sp, time) ;
  }
};

// To prevent multiple definitions, I DECLARE the function in the header only. I then define them
// in the cpp file.
Spatiotemprange sptimeDistance(const Eigen::ArrayXd &, const double &, const Eigen::ArrayXd &,
                               const double &, const std::string &) ;

Eigen::SparseMatrix<double, Eigen::RowMajor> createBlockMatrix(std::vector<Eigen::MatrixXd *>) ;
Eigen::SparseMatrix<double, Eigen::RowMajor> createBlockMatrix(std::vector<Eigen::MatrixXd>) ;
// double logDetBlockMatrix(const Eigen::SparseMatrix<double> &, const Eigen::VectorXi &) ;
double logNormPDF(const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::VectorXd &) ;
double maternCov(const double &, const double &, const double &, const double &, const double &) ;
// Inlining the function solves the linking issue I encountered. As it turns out, templated functions
// must be defined in the header. A simple declaration will lead to a linking error.

template <typename T> inline void deallocate_container(T& c) {
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
};

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> rep(const Eigen::ArrayBase<Derived> & x, const uint times) {
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> container(x.size() * times, 1) ;
  uint index = 0 ;
  for (uint i = 0 ; i < times ; i++) {
    container.segment(index, x.size()) = x ;
    index += x.size() ;
  }
  return container ;
}

template <typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> rep_each(const Eigen::ArrayBase<Derived> & x, const uint times) {
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> container(x.size() * times, 1) ; // double is the most general type. Will
  int index = 0 ;
  for (uint i = 0 ; i < x.size(); i++) {
    container.segment(index, times) = x(i, 0) * Eigen::ArrayBase<Derived>::Ones(times, 1) ;
    index += times ;
  }
  return container ;
}

// Eigen::Matrix<bool, Eigen::Dynamic, 1> operator==(const Eigen::Ref<const Eigen::VectorXi> &, const uint) ;

template<typename Derived>
double median(const Eigen::ArrayBase<Derived> & EigenVec) {
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> VecCopy = EigenVec.matrix() ;
  std::sort(VecCopy.data(), VecCopy.data() + VecCopy.size()) ; // This sorts the vector in descending order, but it doesn't matter for the median!
  double output ;
  int lowerIndex = std::floor(double(VecCopy.size() - 1)/2) ;
  if ((VecCopy.size() % 2) == 1) {
    output = double(VecCopy(lowerIndex)) ;
  } else {
    output = double(VecCopy(lowerIndex) + VecCopy(lowerIndex + 1))/2 ;
  }
  return output ;
}

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
  elem(
    const Eigen::ArrayBase<Derived> & vector,
    const Eigen::Ref<const Eigen::ArrayXi> & indices
  ) {
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> subVector(indices.size(), 1) ;
  for (uint i = 0; i < indices.size(); i++) {
    subVector(i, 1) = vector(indices(i)) ;
  }
  return subVector ;
}

// template<typename Derived>
// Eigen::Matrix<bool, Eigen::Dynamic, 1> operator>(const Eigen::MatrixBase<Derived> & aVector, const uint aConstant) {
//   Eigen::Matrix<bool, Eigen::Dynamic, 1> container(aVector.size()) ;
//   for (uint i = 0; i < aVector.size(); i++) {
//     container(i) = aVector(i) > aConstant ;
//   }
//   return container ;
// }

template<typename Derived>
Eigen::Matrix<bool, Eigen::Dynamic, 1> operator<=(const Eigen::MatrixBase<Derived> & aVector, const uint aConstant) {
  Eigen::Matrix<bool, Eigen::Dynamic, 1> container(aVector.size()) ;
  for (uint i = 0; i < aVector.size(); i++) {
    container(i) = aVector(i) <= aConstant ;
  }
  return container ;
}

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> join_rows(const Eigen::ArrayBase<Derived> & mat1, const Eigen::ArrayBase<Derived> & mat2) {
  if (mat1.rows() != mat2.rows()) {
    throw Rcpp::exception("Error in join_rows: Numbers of cols do not match! \n") ;
  }
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> C(mat1.rows(), mat1.cols() + mat2.cols()) ;
  C << mat1, mat2 ;
  return C ;
}

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> join_cols(const Eigen::ArrayBase<Derived> & mat1, const Eigen::ArrayBase<Derived> & mat2) {
  if (mat1.cols() != mat2.cols()) {
    throw Rcpp::exception("Error in join_cols: Numbers of cols do not match! \n") ;
  }
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> C(mat1.rows() + mat2.rows(), mat1.cols()) ;
  C << mat1,
       mat2 ; // For readability only.
  return C ;
}

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> cols(const Eigen::ArrayBase<Derived> & matrixToSubset, const Eigen::Ref<const Eigen::ArrayXi> & indexVector) {
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> result =  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(matrixToSubset.rows(), indexVector.size()) ;
  for (uint i = 0; i < indexVector.size(); i++) {
    result.col(i) = matrixToSubset.col(indexVector(i)) ;
  }
  return result ;
}

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> rows(const Eigen::ArrayBase<Derived> & matrixToSubset, const Eigen::Ref<const Eigen::ArrayXi> & indexVector) {
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> resultTrans = cols(matrixToSubset.transpose(), indexVector) ;
  return resultTrans.transpose() ;
}

double deg2rad(double) ;
double haversine_distance(double, double, double, double) ;
double vincenty_distance(double, double, double, double) ;

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1>
  sampleWithoutReplacement(
    const Eigen::ArrayBase<Derived> & vecToSampleFrom,
    const int & numElements,
    const bool & sorted,
    std::mt19937_64 & RNG) {
  std::vector<int> aVector ;
  for (uint i = 0; i < vecToSampleFrom.size(); i++) {
    aVector.push_back(i) ;
  }
  std::shuffle(aVector.begin(), aVector.end(), RNG) ;
  aVector.erase(aVector.begin() + numElements, aVector.end()) ;
  if (sorted) {
    std::sort(aVector.begin(), aVector.end()) ;
  }
  Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> subVec(numElements) ;
  for (uint i = 0; i < subVec.size(); i++) {
    subVec(i) = vecToSampleFrom(aVector.at(i)) ;
  }
  return subVec ;
}

Eigen::ArrayXi find(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, 1>>&) ;

template<typename Derived>
Eigen::ArrayXi
  find(const Eigen::ArrayBase<Derived> & vectorToSearch,
       const Eigen::ArrayBase<Derived> & valuesToFind) {
    std::vector<int> indicesVector ;
    for (uint i = 0; i < valuesToFind.size(); i++) {
      for (uint j = 0; j < vectorToSearch.size(); j++) {
        if (i == vectorToSearch(j)) {
          indicesVector.push_back(j) ;
        }
      }
    }
    std::sort(indicesVector.begin(), indicesVector.end()) ;
    Eigen::ArrayXi indicesArray = Eigen::ArrayXi::Map(&indicesVector[0], indicesVector.size()) ;
    return indicesArray ;
}

#endif /* HELPER_H */
