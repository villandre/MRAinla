#include<RcppArmadillo.h>

#ifndef HELPER_H
#define HELPER_H

struct Spatiotemprange{
  double sp ;
  unsigned int time ;

  Spatiotemprange(double & sp, unsigned int & time) : sp(sp), time(time) { } ;
  Spatiotemprange() { } ;
};

// To prevent multiple definitions, I DECLARE the function in the header only. I then define them
// in the cpp file.
Spatiotemprange sptimeDistance(const arma::vec & spCoor1, const unsigned int & time1, const arma::vec & spCoor2,
                               const unsigned int & time2) ;

arma::sp_mat createSparseMatrix(std::vector<arma::mat *>) ;
std::vector<unsigned int> extractBlockIndices(const arma::sp_mat &) ;
arma::mat invertSymmBlockDiag(const arma::sp_mat &, const std::vector<unsigned int> &) ;

// Inlining the function solves the linking issue I encountered. As it turns out, templated functions
// must be defined in the header. A simple declaration will lead to a linking error.
template <typename T> inline void deallocate_container(T& c) {
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
};

#endif /* HELPER_H */
