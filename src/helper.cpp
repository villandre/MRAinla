#include<cmath>
#include<iostream>
#include<cstdlib>

#include "helper.h"

Spatiotemprange sptimeDistance(const arma::vec & spCoor1, const unsigned int & time1, const arma::vec & spCoor2,
                               const unsigned int & time2) {
  arma::vec diffVec = spCoor1 - spCoor2 ;
  arma::vec scaledVec = arma::pow(diffVec, 2) ;
  double sp = arma::sum(scaledVec) ;
  sp = std::sqrt(sp) ;
  int timeDiff = time2 - time1 ;
  timeDiff = abs(timeDiff) ;
  unsigned int timeDiffUint = (unsigned int) timeDiff ;
  return Spatiotemprange(sp, timeDiffUint) ;
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
