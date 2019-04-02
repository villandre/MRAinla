#include<RcppArmadillo.h>

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
Spatiotemprange sptimeDistance(const arma::vec & spCoor1, const double & time1, const arma::vec & spCoor2,
                               const double & time2) ;

arma::sp_mat createBlockMatrix(std::vector<arma::mat *>) ;
arma::uvec extractBlockIndices(const arma::sp_mat &) ;
double logDetBlockMatrix(const arma::sp_mat &) ;
arma::sp_mat invertSymmBlockDiag(const arma::sp_mat &, const arma::uvec &) ;
double logNormPDF(const arma::vec &, const arma::vec &, const arma::vec &) ;
double maternCov(const double &, const double &, const double &, const double &, const double &) ;
// Inlining the function solves the linking issue I encountered. As it turns out, templated functions
// must be defined in the header. A simple declaration will lead to a linking error.
template <typename T> inline void deallocate_container(T& c) {
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
};

#endif /* HELPER_H */
