// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// setupGridCpp
List setupGridCpp(NumericVector responseValues, NumericMatrix spCoords, NumericVector obsTime, NumericMatrix covariateMatrix, uint M, NumericVector lonRange, NumericVector latRange, NumericVector timeRange, uint randomSeed, uint cutForTimeSplit, bool splitTime, int numKnotsRes0, int J);
RcppExport SEXP _MRAinla_setupGridCpp(SEXP responseValuesSEXP, SEXP spCoordsSEXP, SEXP obsTimeSEXP, SEXP covariateMatrixSEXP, SEXP MSEXP, SEXP lonRangeSEXP, SEXP latRangeSEXP, SEXP timeRangeSEXP, SEXP randomSeedSEXP, SEXP cutForTimeSplitSEXP, SEXP splitTimeSEXP, SEXP numKnotsRes0SEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type responseValues(responseValuesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type spCoords(spCoordsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type obsTime(obsTimeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type covariateMatrix(covariateMatrixSEXP);
    Rcpp::traits::input_parameter< uint >::type M(MSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lonRange(lonRangeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type latRange(latRangeSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type timeRange(timeRangeSEXP);
    Rcpp::traits::input_parameter< uint >::type randomSeed(randomSeedSEXP);
    Rcpp::traits::input_parameter< uint >::type cutForTimeSplit(cutForTimeSplitSEXP);
    Rcpp::traits::input_parameter< bool >::type splitTime(splitTimeSEXP);
    Rcpp::traits::input_parameter< int >::type numKnotsRes0(numKnotsRes0SEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(setupGridCpp(responseValues, spCoords, obsTime, covariateMatrix, M, lonRange, latRange, timeRange, randomSeed, cutForTimeSplit, splitTime, numKnotsRes0, J));
    return rcpp_result_gen;
END_RCPP
}
// LogJointHyperMarginalToWrap
double LogJointHyperMarginalToWrap(SEXP treePointer, Rcpp::List MRAhyperparas, double fixedEffSD, double errorSD, Rcpp::List MRAcovParasGammaAlphaBeta, Rcpp::NumericVector FEmuVec, NumericVector fixedEffGammaAlphaBeta, NumericVector errorGammaAlphaBeta, bool matern, double spaceNuggetSD, double timeNuggetSD, bool recordFullConditional, Rcpp::Function optimFun, Rcpp::Function gradCholeskiFun, Rcpp::Function sparseMatrixConstructFun, Rcpp::Function sparseDeterminantFun);
RcppExport SEXP _MRAinla_LogJointHyperMarginalToWrap(SEXP treePointerSEXP, SEXP MRAhyperparasSEXP, SEXP fixedEffSDSEXP, SEXP errorSDSEXP, SEXP MRAcovParasGammaAlphaBetaSEXP, SEXP FEmuVecSEXP, SEXP fixedEffGammaAlphaBetaSEXP, SEXP errorGammaAlphaBetaSEXP, SEXP maternSEXP, SEXP spaceNuggetSDSEXP, SEXP timeNuggetSDSEXP, SEXP recordFullConditionalSEXP, SEXP optimFunSEXP, SEXP gradCholeskiFunSEXP, SEXP sparseMatrixConstructFunSEXP, SEXP sparseDeterminantFunSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type treePointer(treePointerSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MRAhyperparas(MRAhyperparasSEXP);
    Rcpp::traits::input_parameter< double >::type fixedEffSD(fixedEffSDSEXP);
    Rcpp::traits::input_parameter< double >::type errorSD(errorSDSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type MRAcovParasGammaAlphaBeta(MRAcovParasGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type FEmuVec(FEmuVecSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type fixedEffGammaAlphaBeta(fixedEffGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type errorGammaAlphaBeta(errorGammaAlphaBetaSEXP);
    Rcpp::traits::input_parameter< bool >::type matern(maternSEXP);
    Rcpp::traits::input_parameter< double >::type spaceNuggetSD(spaceNuggetSDSEXP);
    Rcpp::traits::input_parameter< double >::type timeNuggetSD(timeNuggetSDSEXP);
    Rcpp::traits::input_parameter< bool >::type recordFullConditional(recordFullConditionalSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type optimFun(optimFunSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type gradCholeskiFun(gradCholeskiFunSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type sparseMatrixConstructFun(sparseMatrixConstructFunSEXP);
    Rcpp::traits::input_parameter< Rcpp::Function >::type sparseDeterminantFun(sparseDeterminantFunSEXP);
    rcpp_result_gen = Rcpp::wrap(LogJointHyperMarginalToWrap(treePointer, MRAhyperparas, fixedEffSD, errorSD, MRAcovParasGammaAlphaBeta, FEmuVec, fixedEffGammaAlphaBeta, errorGammaAlphaBeta, matern, spaceNuggetSD, timeNuggetSD, recordFullConditional, optimFun, gradCholeskiFun, sparseMatrixConstructFun, sparseDeterminantFun));
    return rcpp_result_gen;
END_RCPP
}
// GetFullCondMean
arma::vec GetFullCondMean(SEXP treePointer);
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
arma::vec GetFullCondSDs(SEXP treePointer);
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

static const R_CallMethodDef CallEntries[] = {
    {"_MRAinla_setupGridCpp", (DL_FUNC) &_MRAinla_setupGridCpp, 13},
    {"_MRAinla_LogJointHyperMarginalToWrap", (DL_FUNC) &_MRAinla_LogJointHyperMarginalToWrap, 16},
    {"_MRAinla_GetFullCondMean", (DL_FUNC) &_MRAinla_GetFullCondMean, 1},
    {"_MRAinla_GetFullCondSDs", (DL_FUNC) &_MRAinla_GetFullCondSDs, 1},
    {"_MRAinla_ComputeCondPredStats", (DL_FUNC) &_MRAinla_ComputeCondPredStats, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRAinla(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
