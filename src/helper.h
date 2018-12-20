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

// Inlining the function solves the linking issue I encountered. As it turns out, templated functions
// must be defined in the header. A simple declaration will lead to a linking error.
template <typename T> inline void deallocate_container(T& c) {
  for (typename T::iterator i = c.begin(); i != c.end(); ++i)
    delete *i;
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
    X.submat( idx, idx, idx + dimvec[i] - 1, idx + dimvec[i] - 1 ) = *(listOfMatrices.at(i)) ;
    idx = idx + dimvec[i] ;
  }

  return(X);
}

#endif /* HELPER_H */
