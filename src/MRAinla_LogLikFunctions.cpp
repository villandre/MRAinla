// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <boost/date_time.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <string>
#include <algorithm>
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

namespace bt = boost::posix_time;

const std::locale formats[] = {    // this shows a subset only, see the source file for full list
  std::locale(std::locale::classic(), new bt::time_input_facet("%Y-%m-%d %H:%M:%S%f")),
  std::locale(std::locale::classic(), new bt::time_input_facet("%Y/%m/%d %H:%M:%S%f")),

  std::locale(std::locale::classic(), new bt::time_input_facet("%Y-%m-%d")),
  std::locale(std::locale::classic(), new bt::time_input_facet("%b/%d/%Y")),
};
const size_t nformats = sizeof(formats)/sizeof(formats[0]);

typedef unsigned int uint ;

// [[Rcpp::export]]

SEXP setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, DateVector obsTime, uint M, NumericVector lonRange, NumericVector latRange, DateVector timeRange, uint randomSeed, uint cutForTimeSplit = 200, uint cutForLatSplit = 200)
{
  // TO_DO
}

// [[Rcpp::export]]

List logLikCpp(SEXP treePointer, uint brickDepth)
{
  // TO_DO
}
