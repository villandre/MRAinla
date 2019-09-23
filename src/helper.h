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

Eigen::SparseMatrix<float> createBlockMatrix(std::vector<Eigen::MatrixXf *>) ;
Eigen::VectorXi extractBlockIndices(const Eigen::SparseMatrix<float> &) ;
double logDetBlockMatrix(const Eigen::SparseMatrix<float> &, const Eigen::VectorXi &) ;
Eigen::SparseMatrix<float> invertSymmBlockDiag(const Eigen::SparseMatrix<float> &, const Eigen::VectorXi &) ;
double logNormPDF(const Eigen::VectorXf &, const Eigen::VectorXf &, const Eigen::VectorXf &) ;
double maternCov(const double &, const double &, const double &, const double &, const double &) ;
// Inlining the function solves the linking issue I encountered. As it turns out, templated functions
// must be defined in the header. A simple declaration will lead to a linking error.

template <typename T> inline void deallocate_container(T& c) {
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
};

template <typename T>
Eigen::Matrix<T, Dynamic, 1> rep(const Eigen::Matrix<T, Dynamic, 1> & x, const uint times) {
  Eigen::Matrix<T, Dynamic, 1> container(x.size() * times) ;
  index = 0 ;
  for (uint i = 0 ; i < times ; i++) {
    container.segment(index, index + x.size()) = x ;
    index += x.size() ;
  }
  return container ;
}

template <typename T>
Eigen::Matrix<T, Dynamic, 1> rep_each(Eigen::Matrix<T, Dynamic, 1> & x, const uint times) {
  Eigen::Matrix<T, Dynamic, 1> container(x.size() * times()) ;
  int index = 0 ;
  for (auto & i : x) {
    container.segment(index, times) = x ;
  }
  return container ;
}

#endif /* HELPER_H */
