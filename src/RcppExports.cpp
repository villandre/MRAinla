// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <RcppGSL.h>
#include <Rcpp.h>

using namespace Rcpp;

// setupGridCpp
List setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, NumericMatrix predCoords, NumericVector obsTime, NumericVector predTime, NumericMatrix covariateMatrix, NumericMatrix predCovariateMatrix, uint Mlon, uint Mlat, NumericVector lonRange, NumericVector latRange, uint randomSeed, int numKnotsRes0, double J, String distMethod);
RcppExport SEXP _MRAinla_setupGridCpp(SEXP responseValuesSEXP, SEXP spCoordsSEXP, SEXP predCoordsSEXP, SEXP obsTimeSEXP, SEXP predTimeSEXP, SEXP covariateMatrixSEXP, SEXP predCovariateMatrixSEXP, SEXP MlonSEXP, SEXP MlatSEXP, SEXP lonRangeSEXP, SEXP latRangeSEXP, SEXP randomSeedSEXP, SEXP numKnotsRes0SEXP, SEXP JSEXP, SEXP distMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type responseValues(responseValuesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type spCoords(spCoordsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type predCoords(predCoordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type obsTime(obsTimeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type predTime(predTimeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariateMatrix(covariateMatrixSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type predCovariateMatrix(predCovariateMatrixSEXP);
    Rcpp::traits::input_parameter< uint >::type Mlon(MlonSEXP);
    Rcpp::traits::input_parameter< uint >::type Mlat(MlatSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lonRange(lonRangeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latRange(latRangeSEXP);
    Rcpp::traits::input_parameter< uint >::type randomSeed(randomSeedSEXP);
    Rcpp::traits::input_parameter< int >::type numKnotsRes0(numKnotsRes0SEXP);
    Rcpp::traits::input_parameter< double >::type J(JSEXP);
    Rcpp::traits::input_parameter< String >::type distMethod(distMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(setupGridCpp(responseValues, spCoords, predCoords, obsTime, predTime, covariateMatrix, predCovariateMatrix, Mlon, Mlat, lonRange, latRange, randomSeed, numKnotsRes0, J, distMethod));
    return rcpp_result_gen;
END_RCPP
}
// LogJointHyperMarginalToWrap
double LogJointHyperMarginalToWrap(SEXP treePointer, Rcpp::List MRAhyperparas, double timeCovPara, double fixedEffSD, double errorSD, Rcpp::List MRAcovParasGammaAlphaBeta, Rcpp::NumericVector FEmuVec, NumericVector fixedEffGammaAlphaBeta, NumericVector errorGammaAlphaBeta, NumericVector timeGammaAlphaBeta, double spaceNuggetSD, bool recordFullConditional, bool processPredictions);
RcppExport SEXP _MRAinla_LogJointHyperMarginalToWrap(SEXP treePointerSEXP, SEXP MRAhyperparasSEXP, SEXP timeCovParaSEXP, SEXP fixedEffSDSEXP, SEXP errorSDSEXP, SEXP MRAcovParasGammaAlphaBetaSEXP, SEXP FEmuVecSEXP, SEXP fixedEffGammaAlphaBetaSEXP, SEXP errorGammaAlphaBetaSEXP, SEXP timeGammaAlphaBetaSEXP, SEXP spaceNuggetSDSEXP, SEXP recordFullConditionalSEXP, SEXP processPredictionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MRAhyperparas(MRAhyperparasSEXP);
    Rcpp::traits::input_parameter< double >::type timeCovPara(timeCovParaSEXP);
    Rcpp::traits::input_parameter< double >::type fixedEffSD(fixedEffSDSEXP);
    Rcpp::traits::input_parameter< double >::type errorSD(errorSDSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MRAcovParasGammaAlphaBeta(MRAcovParasGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type FEmuVec(FEmuVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fixedEffGammaAlphaBeta(fixedEffGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type errorGammaAlphaBeta(errorGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type timeGammaAlphaBeta(timeGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< double >::type spaceNuggetSD(spaceNuggetSDSEXP);
    Rcpp::traits::input_parameter< bool >::type recordFullConditional(recordFullConditionalSEXP);
    Rcpp::traits::input_parameter< bool >::type processPredictions(processPredictionsSEXP);
    rcpp_result_gen = Rcpp::wrap(LogJointHyperMarginalToWrap(treePointer, MRAhyperparas, timeCovPara, fixedEffSD, errorSD, MRAcovParasGammaAlphaBeta, FEmuVec, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, timeGammaAlphaBeta, spaceNuggetSD, recordFullConditional, processPredictions));
    return rcpp_result_gen;
END_RCPP
}
// GetFullCondMean
Eigen::VectorXd GetFullCondMean(SEXP treePointer);
RcppExport SEXP _MRAinla_GetFullCondMean(SEXP treePointerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    rcpp_result_gen = Rcpp::wrap(GetFullCondMean(treePointer));
    return rcpp_result_gen;
END_RCPP
}
// GetFullCondSDs
Eigen::VectorXd GetFullCondSDs(SEXP treePointer);
RcppExport SEXP _MRAinla_GetFullCondSDs(SEXP treePointerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    rcpp_result_gen = Rcpp::wrap(GetFullCondSDs(treePointer));
    return rcpp_result_gen;
END_RCPP
}
// ComputeCondPredStats
Rcpp::List ComputeCondPredStats(SEXP treePointer, NumericMatrix spCoordsForPredict, NumericVector timeForPredict, NumericMatrix covariateMatrixForPredict);
RcppExport SEXP _MRAinla_ComputeCondPredStats(SEXP treePointerSEXP, SEXP spCoordsForPredictSEXP, SEXP timeForPredictSEXP, SEXP covariateMatrixForPredictSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type spCoordsForPredict(spCoordsForPredictSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type timeForPredict(timeForPredictSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariateMatrixForPredict(covariateMatrixForPredictSEXP);
    rcpp_result_gen = Rcpp::wrap(ComputeCondPredStats(treePointer, spCoordsForPredict, timeForPredict, covariateMatrixForPredict));
    return rcpp_result_gen;
END_RCPP
}
// GetNumTips
int GetNumTips(SEXP treePointer);
RcppExport SEXP _MRAinla_GetNumTips(SEXP treePointerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    rcpp_result_gen = Rcpp::wrap(GetNumTips(treePointer));
    return rcpp_result_gen;
END_RCPP
}
// GetPredObsOrder
Eigen::VectorXi GetPredObsOrder(SEXP treePointer);
RcppExport SEXP _MRAinla_GetPredObsOrder(SEXP treePointerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    rcpp_result_gen = Rcpp::wrap(GetPredObsOrder(treePointer));
    return rcpp_result_gen;
END_RCPP
}
// GetHmat
Eigen::SparseMatrix<double> GetHmat(SEXP treePointer);
RcppExport SEXP _MRAinla_GetHmat(SEXP treePointerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    rcpp_result_gen = Rcpp::wrap(GetHmat(treePointer));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRAinla_setupGridCpp", (DL_FUNC) &_MRAinla_setupGridCpp, 15},
    {"_MRAinla_LogJointHyperMarginalToWrap", (DL_FUNC) &_MRAinla_LogJointHyperMarginalToWrap, 13},
    {"_MRAinla_GetFullCondMean", (DL_FUNC) &_MRAinla_GetFullCondMean, 1},
    {"_MRAinla_GetFullCondSDs", (DL_FUNC) &_MRAinla_GetFullCondSDs, 1},
    {"_MRAinla_ComputeCondPredStats", (DL_FUNC) &_MRAinla_ComputeCondPredStats, 4},
    {"_MRAinla_GetNumTips", (DL_FUNC) &_MRAinla_GetNumTips, 1},
    {"_MRAinla_GetPredObsOrder", (DL_FUNC) &_MRAinla_GetPredObsOrder, 1},
    {"_MRAinla_GetHmat", (DL_FUNC) &_MRAinla_GetHmat, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRAinla(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
