#include<cmath>

#include "helper.h"

Spatiotemprange sptimeDistance(arma::vec spCoor1, unsigned int time1, arma::vec spCoor2,
                               unsigned int time2) {
  arma::vec diffVec = spCoor1 - spCoor2 ;
  arma::vec scaledVec = arma::pow(diffVec, 2) ;
  double sp = arma::sum(scaledVec) ;
  sp = std::sqrt(sp) ;
  unsigned int timeDiff = time2 - time1 ;
  return Spatiotemprange(sp, timeDiff) ;
};
