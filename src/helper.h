#include<RcppArmadillo.h>

#ifndef HELPER_H
#define HELPER_H

struct Spatiotemprange{
  double sp ;
  unsigned int time ;

  Spatiotemprange(double & sp, unsigned int & time) : sp(sp), time(time) { } ;
  Spatiotemprange() { } ;
};

// To prevent multiple definitions, I DECLARE the functions in the header only. I then define them
// in the cpp file.
Spatiotemprange sptimeDistance(arma::vec spCoor1, unsigned int time1, arma::vec spCoor2,
                               unsigned int time2) ;

// template <typename T> void deallocate_container(T& c) ;

#endif /* HELPER_H */
