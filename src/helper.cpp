#include<math.h>
#include<iostream>
#include<stdlib.h>
#include<gsl/gsl_sf_gamma.h>
#include <boost/math/special_functions/bessel.hpp>

#include "helper.h"

using namespace boost ;
using namespace std ;
using namespace math ;
typedef Eigen::VectorXd vec ;
typedef Eigen::MatrixXd mat ;
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> sp_mat ;
typedef Eigen::VectorXi uvec ;
typedef Eigen::MatrixXi umat ;
typedef Eigen::Triplet<double> Triplet;

// Other choices for method: "vincenty" and "euclidean"
// The assumption is that spCoor1 is a vector of length 2 which lists longitude first, then latitude.
Spatiotemprange sptimeDistance(const Eigen::ArrayXd & spCoor1, const double & time1, const Eigen::ArrayXd & spCoor2,
                               const double & time2, const string & method) {
  double sp = 0 ;
  if (method == "haversine") {
    sp = haversine_distance(spCoor1(0), spCoor1(1), spCoor2(0), spCoor2(1)) ; // Works under the assumption that coordinates are (longitude, latitude)
  } else { // Input other
    Rcpp::stop("For now, only Haversine distance is implemented. \n") ;
    // sp = vincenty_distance(spCoor1(0), spCoor1(1), spCoor2(0), spCoor2(1)) ; // Same as for haversine
  }
  double timeDiff = abs(time2 - time1) ;

  return Spatiotemprange(sp, timeDiff) ;
};

sp_mat createBlockMatrix(std::vector<mat *> listOfMatrices) {

  std::vector<Triplet> tripletList;

  int offset = 0 ;

  for (auto & aMatrix : listOfMatrices) {
    for (uint i = 0; i < aMatrix->rows() ; i++) {
      for (uint j = 0; j < aMatrix->cols(); j++) {
        tripletList.push_back(Triplet(i + offset, j + offset, (*aMatrix)(i, j))) ;
      }
    }
    offset += aMatrix->rows() ;
  }
  sp_mat X(offset, offset);

  X.setFromTriplets(tripletList.begin(), tripletList.end()) ;

  return X ;
}

sp_mat createBlockMatrix(std::vector<mat> listOfMatrices) {

  std::vector<Triplet> tripletList;

  int offset = 0 ;

  for (auto & aMatrix : listOfMatrices) {
    for (uint i = 0; i < aMatrix.rows() ; i++) {
      for (uint j = 0; j < aMatrix.cols(); j++) {
        tripletList.push_back(Triplet(i + offset, j + offset, aMatrix(i, j))) ;
      }
    }
    offset += aMatrix.rows() ;
  }
  sp_mat X(offset, offset);

  X.setFromTriplets(tripletList.begin(), tripletList.end()) ;

  return X ;
}

double logNormPDF(const vec & x, const vec & mu, const vec & sd) {
  double logValue = 0;
  for (unsigned int i = 0 ; i < x.size() ; i++) {
    logValue += (-log(sd(i)) -
      0.5 * pow((x(i) - mu(i)), 2)/pow(sd(i), 2)) ;
  }
  logValue += (-0.5*x.size()*(log(2) + log(PI))) ;
  return logValue ;
}

// The notation is borrowed from the Wikipedia entry.

double maternCov(const double & distance, const double & rho,
                    const double & smoothness, const double & scale, const double & nugget) {
  double maternValue ;

  if (distance == 0) {
    maternValue = pow(scale, 2) + nugget ;
  } else {
    if (smoothness >= 1e6) {
      maternValue = pow(scale, 2) * exp(-pow(distance, 2)/(2 * pow(rho, 2))) ;
    } else if (smoothness == 0.5) {
      maternValue = pow(scale, 2) * exp(-distance / rho) ;
    } else if (smoothness == 1.5) {
      maternValue = pow(scale, 2) *
        (1 + sqrt(3) * distance / rho) *
        exp(-sqrt(3) * distance / rho) ;
    } else if (smoothness == 2.5) {
      maternValue = pow(scale, 2) *
        (1 + sqrt(5) * distance/rho + 5 * pow(distance, 2)/(3 * pow(rho, 2))) *
        exp(-sqrt(5) * distance/rho) ;
    } else {
      double base = pow(2 * smoothness, 0.5) * distance / rho ;
      double bessel = boost::math::cyl_bessel_k(smoothness, base) ;
      maternValue = pow(scale, 2) * pow(2, 1 - smoothness) / gsl_sf_gamma(smoothness) *
        pow(base, smoothness) * bessel ;
    }
  }
  return maternValue ;
}

Eigen::ArrayXi find(const Eigen::Ref<const Eigen::Array<bool, Eigen::Dynamic, 1>>& logicalVector) {
  Eigen::ArrayXi outputVec(logicalVector.size()) ;
  uint index = 0 ;
  for (uint i = 0 ; i < logicalVector.size(); i++) {
    if (logicalVector(i)) {
      outputVec(index) = i ;
      index += 1 ;
    }
  }
  return outputVec.segment(0, index) ;
}

// The following code was taken from https://gist.github.com/ed-flanagan/e6dc6b8d3383ef5a354a.

/*
 * Great-circle distance computational forumlas
 *
 * https://en.wikipedia.org/wiki/Great-circle_distance
 */

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef M_PI
#define M_PI    3.1415926535897932384626433832795
#endif

using namespace std;

static const double earth_radius_km = 6371.0;


// Obtained from https://www.geeksforgeeks.org/haversine-formula-to-find-distance-between-two-points-on-a-sphere/

double haversine_distance(double lon1, double lat1,
                        double lon2, double lat2)
{
  // distance between latitudes
  // and longitudes
  double dLat = (lat2 - lat1) *
    M_PI / 180.0;
  double dLon = (lon2 - lon1) *
    M_PI / 180.0;

  // convert to radians
  lat1 = (lat1) * M_PI / 180.0;
  lat2 = (lat2) * M_PI / 180.0;

  // apply formulae
  double a = pow(sin(dLat / 2), 2) +
    pow(sin(dLon / 2), 2) *
    cos(lat1) * cos(lat2);

  double c = 2 * asin(sqrt(a));
  return earth_radius_km * c;
}

double vincenty_distance(double longitude1, double latitude1, double longitude2, double latitude2)
{
  Rcpp::stop("Vincenty distance not implemented for now... \n") ;
}


