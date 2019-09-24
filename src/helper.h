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
Spatiotemprange sptimeDistance(const Eigen::VectorXf & spCoor1, const double & time1, const Eigen::VectorXf & spCoor2,
                               const double & time2) ;

Eigen::SparseMatrix<double> createBlockMatrix(std::vector<Eigen::MatrixXd *>) ;
Eigen::VectorXi extractBlockIndices(const Eigen::SparseMatrix<double> &) ;
double logDetBlockMatrix(const Eigen::SparseMatrix<double> &, const Eigen::VectorXi &) ;
Eigen::SparseMatrix<double> invertSymmBlockDiag(const Eigen::SparseMatrix<double> &, const Eigen::VectorXi &) ;
double logNormPDF(const Eigen::VectorXf &, const Eigen::VectorXf &, const Eigen::VectorXf &) ;
double maternCov(const double &, const double &, const double &, const double &, const double &) ;
// Inlining the function solves the linking issue I encountered. As it turns out, templated functions
// must be defined in the header. A simple declaration will lead to a linking error.

template <typename T> inline void deallocate_container(T& c) {
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
};

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> rep(const Eigen::Matrix<T, Eigen::Dynamic, 1> & x, const uint times) {
  Eigen::Matrix<T, Eigen::Dynamic, 1> container(x.size() * times) ;
  index = 0 ;
  for (uint i = 0 ; i < times ; i++) {
    container.segment(index, index + x.size()) = x ;
    index += x.size() ;
  }
  return container ;
}

template <typename Derived>
Eigen::MatrixBase<Derived> rep_each(const Eigen::MatrixBase<Derived> & x, const uint times) {
  Eigen::MatrixBase<Derived> container(x.size() * times, 1) ; // double is the most general type. Will
  int index = 0 ;
  for (auto & i : x) {
    container.segment(index, times) = x ;
  }
  return container ;
}

Eigen::Matrix<bool, Eigen::Dynamic, 1> operator==(const Eigen::Ref<const Eigen::VectorXi> &, const uint) ;

template<typename Derived>
double median(const Eigen::MatrixBase<Derived> & EigenVec) {
  Eigen::MatrixBase<Derived> VecCopy = EigenVec ;
  std::sort(VecCopy.data(), VecCopy.data() + VecCopy.size()) ; // This sorts the vector is descending order, but it doesn't matter for the median!
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
Eigen::MatrixBase<Derived> elem(const Eigen::MatrixBase<Derived> & vector, const Eigen::Ref<const Eigen::VectorXi> & indices) {
  Eigen::MatrixBase<Derived> subVector(indices.size()) ;
  for (uint i = 0; i < indices.size(); i++) {
    subVector(i) = vector(indices(i)) ;
  }
  return subVector ;
}

template<typename Derived>
Eigen::Matrix<bool, Eigen::Dynamic, 1> operator>(const Eigen::MatrixBase<Derived> & aVector, const uint aConstant) {
  Eigen::Matrix<bool, Eigen::Dynamic, 1> container(aVector.size()) ;
  for (uint i = 0; i < aVector.size(); i++) {
    container(i) = aVector(i) > aConstant ;
  }
  return container ;
}

Eigen::Matrix<bool, Eigen::Dynamic, 1> find(const Eigen::Ref<const Eigen::Matrix<bool, Eigen::Dynamic, 1>> &) ;

#endif /* HELPER_H */
