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
typedef Eigen::SparseMatrix<double> sp_mat ;
typedef Eigen::VectorXi uvec ;
typedef Eigen::MatrixXi umat ;
typedef Eigen::Triplet<double> Triplet;

// Other choices for method: "vincenty" and "euclidean"
// The assumption is that spCoor1 is a vector of length 2 which lists latitude first, then longitude.
Spatiotemprange sptimeDistance(const Eigen::ArrayXd & spCoor1, const double & time1, const Eigen::ArrayXd & spCoor2,
                               const double & time2, const string & method="haversine") {
  double sp = 0 ;
  if (method == "euclidean") {
    Eigen::ArrayXd diffVec = spCoor1 - spCoor2 ;
    Eigen::ArrayXd scaledVec = diffVec.pow(2) ;
    double sp = scaledVec.sum() ;
    sp = std::sqrt(sp) ;
  } else if (method == "haversine") {
    sp = haversine_distance(spCoor1(0), spCoor1(1), spCoor2(0), spCoor2(1)) ;
  } else if (method == "vincenty") {
    sp = vincenty_distance(spCoor1(0), spCoor1(1), spCoor2(0), spCoor2(1)) ;
  }
  double timeDiff = abs(time2 - time1) ;

  return Spatiotemprange(sp, timeDiff) ;
};

// Pretty slow. Should not be called too often.

sp_mat createBlockMatrix(std::vector<mat *> listOfMatrices) {
  uint numRows = 0 ;
  uint numCols = 0 ;

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

double deg2rad(double deg)
{
  return (deg * M_PI / 180.0);
}

double haversine_distance(double latitude1, double longitude1, double latitude2,
                          double longitude2)
{
  double lat1 = deg2rad(latitude1);
  double lon1 = deg2rad(longitude1);
  double lat2 = deg2rad(latitude2);
  double lon2 = deg2rad(longitude2);

  double d_lat = abs(lat1 - lat2);
  double d_lon = abs(lon1 - lon2);

  double a = pow(sin(d_lat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(d_lon / 2), 2);

  //double d_sigma = 2 * atan2(sqrt(a), sqrt(1 - a));
  double d_sigma = 2 * asin(sqrt(a));

  return earth_radius_km * d_sigma;
}

double vincenty_distance(double latitude1, double longitude1, double latitude2,
                         double longitude2)
{
  double lat1 = deg2rad(latitude1);
  double lon1 = deg2rad(longitude1);
  double lat2 = deg2rad(latitude2);
  double lon2 = deg2rad(longitude2);

  double d_lon = abs(lon1 - lon2);

  // Numerator
  double a = pow(cos(lat2) * sin(d_lon), 2);

  double b = cos(lat1) * sin(lat2);
  double c = sin(lat1) * cos(lat2) * cos(d_lon);
  double d = pow(b - c, 2);

  double e = sqrt(a + d);

  // Denominator
  double f = sin(lat1) * sin(lat2);
  double g = cos(lat1) * cos(lat2) * cos(d_lon);

  double h = f + g;

  double d_sigma = atan2(e, h);

  return earth_radius_km * d_sigma;
}


