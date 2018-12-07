// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// setupGridCpp
SEXP setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, IntegerVector obsTime, NumericMatrix covariateMatrix, uint M, NumericVector lonRange, NumericVector latRange, IntegerVector timeRange, uint randomSeed, uint cutForTimeSplit);
RcppExport SEXP _MRAinla_setupGridCpp(SEXP responseValuesSEXP, SEXP spCoordsSEXP, SEXP obsTimeSEXP, SEXP covariateMatrixSEXP, SEXP MSEXP, SEXP lonRangeSEXP, SEXP latRangeSEXP, SEXP timeRangeSEXP, SEXP randomSeedSEXP, SEXP cutForTimeSplitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type responseValues(responseValuesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type spCoords(spCoordsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type obsTime(obsTimeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariateMatrix(covariateMatrixSEXP);
    Rcpp::traits::input_parameter< uint >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lonRange(lonRangeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latRange(latRangeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type timeRange(timeRangeSEXP);
    Rcpp::traits::input_parameter< uint >::type randomSeed(randomSeedSEXP);
    Rcpp::traits::input_parameter< uint >::type cutForTimeSplit(cutForTimeSplitSEXP);
    rcpp_result_gen = Rcpp::wrap(setupGridCpp(responseValues, spCoords, obsTime, covariateMatrix, M, lonRange, latRange, timeRange, randomSeed, cutForTimeSplit));
    return rcpp_result_gen;
END_RCPP
}
// logLikMRAcpp
List logLikMRAcpp(SEXP treePointer, NumericVector& covParameters, NumericVector& fixedEffectParameters, double& errorSD, double& fixedEffSD);
RcppExport SEXP _MRAinla_logLikMRAcpp(SEXP treePointerSEXP, SEXP covParametersSEXP, SEXP fixedEffectParametersSEXP, SEXP errorSDSEXP, SEXP fixedEffSDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type covParameters(covParametersSEXP);
    Rcpp::traits::input_parameter< NumericVector& >::type fixedEffectParameters(fixedEffectParametersSEXP);
    Rcpp::traits::input_parameter< double& >::type errorSD(errorSDSEXP);
    Rcpp::traits::input_parameter< double& >::type fixedEffSD(fixedEffSDSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikMRAcpp(treePointer, covParameters, fixedEffectParameters, errorSD, fixedEffSD));
    return rcpp_result_gen;
END_RCPP
}
// predictMRArcpp
List predictMRArcpp(SEXP treePointer, NumericMatrix predSpatialCoor, IntegerVector predTime);
RcppExport SEXP _MRAinla_predictMRArcpp(SEXP treePointerSEXP, SEXP predSpatialCoorSEXP, SEXP predTimeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type predSpatialCoor(predSpatialCoorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type predTime(predTimeSEXP);
    rcpp_result_gen = Rcpp::wrap(predictMRArcpp(treePointer, predSpatialCoor, predTime));
    return rcpp_result_gen;
END_RCPP
}
// inla
NumericMatrix inla(SEXP treePointer, NumericMatrix predictionLocations, IntegerVector predictionTime, double stepSize);
RcppExport SEXP _MRAinla_inla(SEXP treePointerSEXP, SEXP predictionLocationsSEXP, SEXP predictionTimeSEXP, SEXP stepSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type predictionLocations(predictionLocationsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type predictionTime(predictionTimeSEXP);
    Rcpp::traits::input_parameter< double >::type stepSize(stepSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(inla(treePointer, predictionLocations, predictionTime, stepSize));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRAinla_setupGridCpp", (DL_FUNC) &_MRAinla_setupGridCpp, 10},
    {"_MRAinla_logLikMRAcpp", (DL_FUNC) &_MRAinla_logLikMRAcpp, 5},
    {"_MRAinla_predictMRArcpp", (DL_FUNC) &_MRAinla_predictMRArcpp, 3},
    {"_MRAinla_inla", (DL_FUNC) &_MRAinla_inla, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRAinla(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
