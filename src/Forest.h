#include "AugTree.h"

#ifndef MYPROJECT_FOREST_H
#define MYPROJECT_FOREST_H

namespace MRAinla {
class Forest {
public:
  Forest(uint &, uint &, Eigen::Array2d &, Eigen::Array2d &, vec &, Eigen::ArrayXXd &, Eigen::ArrayXd &, Eigen::ArrayXXd &, Eigen::ArrayXXd &, Eigen::ArrayXd &, unsigned long int &, Eigen::ArrayXXd &, const unsigned int, double, const std::string &) ;

  void SetFixedEffParameters(vec & fixedParas) {
    if (fixedParas.isApprox(m_fixedEffParameters)) m_recomputeGlobalLogLik = true ;
    m_fixedEffParameters = fixedParas ;
  }
  void SetFixedEffSD(const double & fixedEffSD) {
    m_fixedEffSD = fixedEffSD ;
  }
  void SetErrorSD(const double & errorSD) {
    if (errorSD != m_errorSD) {
      m_recomputeGlobalLogLik = true ;
    }
    m_errorSD = errorSD ;
  }
  void SetMRAcovParas(const Rcpp::List &) ;
  void SetProcessPredictions(bool processPredictions) {
    m_processPredictions = processPredictions ;
  }
  void SetErrorGammaAlphaBeta(const GammaHyperParas & alphaBeta) {m_errorGammaAlphaBeta = alphaBeta ;}
  void SetFixedEffGammaAlphaBeta(const GammaHyperParas & alphaBeta) {m_fixedEffGammaAlphaBeta = alphaBeta ;}
  void SetTimeCovParaGammaAlphaBeta(const GammaHyperParas & alphaBeta) {m_timeCovParaGammaAlphaBeta = alphaBeta ;}
  void SetMRAcovParasGammaAlphaBeta(const Rcpp::List &) ;
  void SetFEmu(const vec & muVec) {m_FEmu = muVec ;}
  void SetSpaceNuggetSD(const double & spaceNuggetSD) {
    m_spaceNuggetSD = spaceNuggetSD ;
  }
  void SetRecordFullConditional(const bool recordIt) { m_recordFullConditional = recordIt ;}
  void SetPredictData(const mat & spCoords, const vec & timeValues, const mat covariates) {
    vec placeholder(covariates.rows()) ;
    inputdata dataObject(placeholder, spCoords, timeValues, covariates) ;
    m_predictData = dataObject ;
  }
  void ToggleGammaParasSet() { m_GammaParasSet = true ;}

  bool CheckMRAcovParasGammaAlphaBeta() {
    return m_GammaParasSet ;
  }
  void ComputeLogJointPsiMarginal() ;
  inputdata GetDataset() {return m_dataset;}
  double GetLogJointPsiMarginal() { return m_logJointPsiMarginal ;}
  vec GetFullCondMean() { return m_FullCondMean ;}
  vec GetFullCondSDs() { return m_FullCondSDs ;}
  sp_mat & GetHmatPred() {return m_HmatPred ;}
  sp_mat & GetHmat() {return m_Hmat ;}
  Eigen::ArrayXi GetObsOrderForFpredMat() { return m_obsOrderForFpredMat ;}

  vec ComputeEvar() ;

private:

  std::vector<AugTree> m_treeVector ;
  bool m_recomputeMRAlogLik{ true } ; // When this flag is true, more computations are required to get the log-lik.
  bool m_recomputeGlobalLogLik{ true } ; // When this flag is true, the global log-likelihood (conditional on all mean parameters) needs to be recomputed.
  bool m_processPredictions{ false } ; // When this flag is true, the H matrix for predictions is computed.

  double m_globalLogLik{ 0 } ;
  double m_logPrior{ 0 } ;
  double m_logCondDist{ 0 } ;
  double m_logFullCond{ 0 } ;
  double m_logJointPsiMarginal{ 0 } ;

  double m_spaceNuggetSD{ 1e-6 } ;
  std::string m_distMethod ;
  double m_errorSD{ 0 } ;
  double m_fixedEffSD{ 0 } ;
  maternGammaPriorParasWithoutScale m_maternParasGammaAlphaBetaSpace ;
  GammaHyperParas m_maternSpacetimeScalingGammaAlphaBeta ;
  GammaHyperParas m_fixedEffGammaAlphaBeta ;
  GammaHyperParas m_errorGammaAlphaBeta ;
  GammaHyperParas m_timeCovParaGammaAlphaBeta ;
  maternVec m_MRAcovParasSpace ;
  spaceDimensions m_mapDimensions;

  inputdata m_dataset ; // First element is response, second is spatial coordinates, last is time.

  std::mt19937_64 m_randomNumGenerator ;

  double m_spacetimeScaling{ 1e-5 } ;
  vec m_fixedEffParameters ;
  vec m_FEmu ;
  double m_timeCovPara{ 0 } ;
  vec m_Vstar ;
  vec m_FullCondMean ;
  vec m_FullCondSDs ;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_FullCondPrecisionChol ;

  bool m_recordFullConditional{ false } ;
  bool m_GammaParasSet{ false } ;

  inputdata m_predictData ;
  Eigen::ArrayXi m_obsOrderForFpredMat ;

  sp_mat m_Hmat ;
  Eigen::ArrayXi m_obsOrderForFmat ;
  sp_mat m_SigmaFEandEtaInv ;
  Eigen::ArrayXi m_obsOrderForHmat ;
  std::vector<pointerOffset> m_pointerOffsetForHmat ;
  std::vector<pointerOffset> m_pointerOffsetForHmatPred ;

  sp_mat m_HmatPred ;
  sp_mat m_SigmaBetaEtaInv ;
  Eigen::ArrayXd m_uniqueTimeValues ;
  Eigen::Array<bool, Eigen::Dynamic, 1> m_assignedPredToKnot ;

  uint m_numKnots{ 0 } ;

  vec ComputeFullConditionalMean(const vec &) ;
  void ComputeFullCondSDsFE() ;

  std::vector<GaussDistParas> ComputeConditionalPrediction(const inputdata &) ;
  void ComputeLogFCandLogCDandDataLL() ;
  void ComputeLogPriors() ;

  void createHmatrixPred() ;
  void CreateSigmaBetaEtaInvMat() ;
  void UpdateSigmaBetaEtaInvMat() ;

  void ComputeHpred() ;
  void createHmatrix() ;
  void updateHmatrix() ;
  void updateHmatrixPred() ;
  mat createTimePrecisionMatrix(const double, const Eigen::ArrayXd) { };
};
}
#endif
