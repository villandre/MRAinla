#include<RcppEigen.h>

#ifndef HELPER_H
#define HELPER_H

struct Spatiotemprange{
  double sp ;
  double time ;

  Spatiotemprange(double & sp, double & time) : sp(sp), time(time) { } ;
  Spatiotemprange() { } ;
};

// To prevent multiple definitions, I DECLARE the function in the header only. I then define them
// in the cpp file.
Spatiotemprange sptimeDistance(const Eigen::VectorXd & spCoor1, const double & time1, const Eigen::VectorXd & spCoor2,
                               const double & time2) ;

Eigen::SparseMatrix<double> createBlockMatrix(std::vector<Eigen::MatrixXd *>) ;

// double logDetBlockMatrix(const Eigen::SparseMatrix<double> &, const Eigen::VectorXi &) ;
Eigen::SparseMatrix<double> invertSymmBlockDiag(const Eigen::SparseMatrix<double> &, const Eigen::VectorXi &) ;
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
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, 1> VecCopy = EigenVec ;
  std::sort(VecCopy.data(), VecCopy.data() + VecCopy.size()) ; // This sorts the vector in descending order, but it doesn't matter for the median!
  double output ;
  int lowerIndex = std::floor(double(VecCopy.size())/2) ;
  if ((VecCopy.size() % 2) == 1) {
    output = double(VecCopy(lowerIndex)) ;
  } else {
    output = double(VecCopy(lowerIndex) + VecCopy(lowerIndex + 1))/2 ;
  }
  return output ;
}

template<typename Derived>
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> elem(const Eigen::ArrayBase<Derived> & vector, const Eigen::Ref<const Eigen::ArrayXi> & indices) {
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

Eigen::ArrayXi find(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, 1>> &) ;

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

#endif /* HELPER_H */
