#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

struct GammaHyperParas{
  // Default IG mean is 1 (mean = beta/(alpha - 1)), var is beta^2/[(alpha-1)^2(alpha-2)]
  // If alpha <= 2, variance is infinite,
  double m_alpha = 3;
  double m_beta = 2;
  GammaHyperParas() {}
  GammaHyperParas(double alpha, double beta) : m_alpha(alpha), m_beta(beta) { }
  GammaHyperParas(const Eigen::VectorXd & alphaBeta) {
    m_alpha = alphaBeta(0) ;
    m_beta = alphaBeta(1) ;
  }
};

struct maternGammaPriorParasWithoutScale{
  GammaHyperParas m_rho ;
  GammaHyperParas m_smoothness ;

  maternGammaPriorParasWithoutScale() { }
  maternGammaPriorParasWithoutScale(const GammaHyperParas & rho, const GammaHyperParas & smoothness) : m_rho(rho), m_smoothness(smoothness) { }
};

struct maternGammaPriorParas : public maternGammaPriorParasWithoutScale{
  GammaHyperParas m_scale ;

  maternGammaPriorParas(const GammaHyperParas & rho, const GammaHyperParas & smoothness, const GammaHyperParas & scale) : maternGammaPriorParasWithoutScale(rho, smoothness), m_scale(scale) { } ;
};


namespace MRAinla {

class AugTree
{
public:
  AugTree(uint &, uint &, uint &, Eigen::Array2d &, Eigen::Array2d &, Eigen::Array2d &, vec &, Eigen::ArrayXXd &, Eigen::ArrayXd &, Eigen::ArrayXXd &, Eigen::ArrayXXd &, Eigen::ArrayXd &, uint &, unsigned long int &, Eigen::ArrayXXd &, const bool, const unsigned int, double, const std::string &) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

  std::vector<GaussDistParas> ComputeConditionalPrediction(const inputdata &) ;

  void ComputeLogFCandLogCDandDataLL() ;
  void ComputeLogPriors() ;

  // double GetMRAlogLik() const {return m_MRAlogLik ;}
  inputdata GetDataset() {return m_dataset;}
  int GetNumTips() {return m_numTips ;}
  int GetNumKnots() {return m_numKnots ;}
  sp_mat & GetHmatPred() {return m_HmatPred ;}
  sp_mat & GetHmat() {return m_Hmat ;}

  int GetM() { return m_M ;}
  double GetLogJointPsiMarginal() { return m_logJointPsiMarginal ;}
  Eigen::ArrayXi GetObsOrderForHpredMat() { return m_obsOrderForHpredMat ;}

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

  std::vector<TreeNode *> GetLevelNodes(const uint & level) ;
  std::vector<TreeNode *> GetTipNodes() { return GetLevelNodes(m_M) ;}

  void SetErrorGammaAlphaBeta(const GammaHyperParas & alphaBeta) {m_errorGammaAlphaBeta = alphaBeta ;}
  void SetFixedEffGammaAlphaBeta(const GammaHyperParas & alphaBeta) {m_fixedEffGammaAlphaBeta = alphaBeta ;}
  void SetMRAcovParasGammaAlphaBeta(const Rcpp::List &) ;
  void SetFEmu(const vec & muVec) {m_FEmu = muVec ;}
  void SetSpaceAndTimeNuggetSD(const double & spaceNuggetSD, const double & timeNuggetSD) {
    m_spaceNuggetSD = spaceNuggetSD ;
    m_timeNuggetSD = timeNuggetSD ;
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

  void CleanPredictionComponents() ;

  void createHmatrixPred() ;
  void CreateSigmaBetaEtaInvMat() ;
  void UpdateSigmaBetaEtaInvMat() ;
  // sp_mat createQ() ;

  void ComputeLogJointPsiMarginal() ;
  vec GetFullCondMean() { return m_FullCondMean ;}
  vec GetFullCondSDs() { return m_FullCondSDs ;}

  void ComputeHpred() ;
  vec ComputeEvar() ;

  ~ AugTree() {
    deallocate_container(m_vertexVector) ;}

private:

  std::vector<TreeNode *> m_vertexVector ;

  int GetNodePos(int nodeId) {
    int nodePos = 0 ;
    while (m_vertexVector.at(nodePos)->GetNodeId() != nodeId) {
      nodePos += 1 ;
    }
  return nodePos ;
  }

  bool m_recomputeMRAlogLik{ true } ; // When this flag is true, more computations are required to get the log-lik.
  bool m_recomputeGlobalLogLik{ true } ; // When this flag is true, the global log-likelihood (conditional on all mean parameters) needs to be recomputed.
  bool m_processPredictions{ false } ; // When this flag is true, the H matrix for predictions is computed.

  double m_globalLogLik{ 0 } ;
  double m_logPrior{ 0 } ;
  double m_logCondDist{ 0 } ;
  double m_logFullCond{ 0 } ;
  double m_logJointPsiMarginal{ 0 } ;
  int m_numKnots{ 0 } ;
  double m_spaceNuggetSD{ 1e-6 };
  double m_timeNuggetSD{ 1e-6 } ;
  std::string m_distMethod ;

  int m_M{ 0 } ;
  int m_Mlon{ 0 } ;
  int m_Mlat{ 0 } ;
  int m_Mtime{ 0 } ;
  int m_numTips{ 0 } ;
  double m_errorSD{ 0 } ;
  double m_fixedEffSD{ 0 } ;

  maternGammaPriorParasWithoutScale m_maternParasGammaAlphaBetaSpace ;
  maternGammaPriorParasWithoutScale m_maternParasGammaAlphaBetaTime ;
  GammaHyperParas m_maternSpacetimeScalingGammaAlphaBeta ;
  GammaHyperParas m_fixedEffGammaAlphaBeta ;
  GammaHyperParas m_errorGammaAlphaBeta ;

  maternVec m_MRAcovParasSpace ;
  maternVec m_MRAcovParasTime ;
  double m_spacetimeScaling{ 1e-5 } ;
  vec m_fixedEffParameters ;
  vec m_FEmu ;

  dimensions m_mapDimensions;

  inputdata m_dataset ; // First element is response, second is spatial coordinates, last is time.

  std::mt19937_64 m_randomNumGenerator ;

  std::vector<TreeNode *> Descendants(std::vector<TreeNode *>) ;
  void diveAndUpdate(TreeNode *, std::vector<TreeNode *> *) ;

  // Tree construction functions //
  void BuildTree(const uint &, const bool, const unsigned int, double) ;
  void createLevels(TreeNode *, std::string, Eigen::ArrayXi) ;
  void generateKnots(TreeNode *, const unsigned int, double) ;

  void numberNodes() ;
  void computeWmats() ;
  std::vector<TreeNode *> GetLevel(const uint) ;
  Eigen::ArrayXi m_obsOrderForFmat ;

  Eigen::Array<bool, Eigen::Dynamic, 1> m_assignedPredToKnot ;

  inputdata m_predictData ;
  Eigen::ArrayXi m_obsOrderForHpredMat ;
  vec m_Vstar ;
  vec m_FullCondMean ;
  vec m_FullCondSDs ;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> m_FullCondPrecisionChol ;

  sp_mat m_Hmat ;
  std::vector<pointerOffset> m_pointerOffsetForHmat ;
  std::vector<pointerOffset> m_pointerOffsetForHmatPred ;

  sp_mat m_HmatPred ;
  sp_mat m_SigmaFEandEtaInv ;
  bool m_recordFullConditional{ false } ;
  bool m_GammaParasSet{ false } ;

  std::vector<mat *> getKmatricesInversePointers() {
    std::vector<mat *> KmatrixInverseList(m_vertexVector.size()) ;
    for (auto & i : m_vertexVector) KmatrixInverseList.at(i->GetNodeId()) = i->GetKmatrixInverseAddress() ;
    return KmatrixInverseList ;
  }

  void createHmatrix() ;
  void updateHmatrix() ;
  void updateHmatrixPred() ;
  vec ComputeFullConditionalMean(const vec &) ;
  void ComputeFullCondSDsFE() ;
  template <typename MatrixType>
  inline typename MatrixType::Scalar logdet(const MatrixType& M, bool use_cholesky = false) ;

};
}
#endif
