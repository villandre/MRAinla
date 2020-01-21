#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

namespace MRAinla {
struct TwoParsProbDist{
  double m_firstPar ;
  double m_secondPar ;

  void baseInitialise(const Eigen::ArrayXd & twoParsArray) {
    m_firstPar = twoParsArray(0) ;
    m_secondPar = twoParsArray(1) ;
  }
  void baseInitialise(const double & a, const double & b) {
    m_firstPar = a ;
    m_secondPar = b ;
  }
  virtual double computeLogDensity(const double x) const = 0;

  virtual ~ TwoParsProbDist() { } ;
};

struct GammaDist : public TwoParsProbDist {
  GammaDist() { }
  GammaDist(double alpha, double beta) {
    baseInitialise(alpha, beta) ;
  }

  GammaDist(const Eigen::ArrayXd & alphaBeta) {
    baseInitialise(alphaBeta) ;
  }

  double computeLogDensity(const double x)  const override {
    double value = (m_firstPar - 1) * log(x) - x * m_secondPar + m_firstPar * log(m_secondPar) - lgamma(m_firstPar) ;
    return value ;
  }
};

struct NormalDist : public TwoParsProbDist {
  NormalDist() { }
  NormalDist(double mu, double sigma) {
    baseInitialise(mu, sigma) ;
  }
  NormalDist(const Eigen::ArrayXd & muSigma) {
    baseInitialise(muSigma) ;
  }

  double computeLogDensity(const double x) const override {
    double value = -log(m_secondPar) - 0.5 * (log(2) + log(PI)) - 0.5 * pow(x - m_firstPar, 2) * pow(m_secondPar, -2) ;
    return value ;
  }
};

class AugTree
{
public:
  AugTree(const uint &,
          const uint &,
          const uint &,
          const Eigen::Array2d &,
          const Eigen::Array2d &,
          const Eigen::Array2d &,
          const vec &,
          const Eigen::ArrayXXd &,
          const Eigen::ArrayXd &,
          const Eigen::ArrayXXd &,
          const Eigen::ArrayXXd &,
          const Eigen::ArrayXd &,
          const uint &,
          const unsigned long int &,
          const Eigen::ArrayXXd &,
          const bool &,
          const unsigned int &,
          const double &,
          const std::string &,
          const Rcpp::List &,
          const Rcpp::NumericVector &,
          const Rcpp::NumericVector &,
          const Rcpp::NumericVector &,
          const double &,
          const bool &) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

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
  void SetMaternPars(const Rcpp::List &) ;

  void SetProcessPredictions(bool processPredictions) {
    m_processPredictions = processPredictions ;
  }

  std::vector<TreeNode *> GetLevelNodes(const uint & level) ;
  std::vector<TreeNode *> GetTipNodes() { return GetLevelNodes(m_M) ;}

  void SetParsHyperpars(const Rcpp::List &, const Rcpp::NumericVector &, const Rcpp::NumericVector &, const bool) ;

  void SetFEmu(const vec & muVec) {m_FEmu = muVec ;}
  void SetNuggetSD(const double & nuggetSD) {
    m_nuggetSD = nuggetSD ;
  }
  void SetRecordFullConditional(const bool recordIt) { m_recordFullConditional = recordIt ;}

  void SetPredictData(const mat & spCoords, const vec & timeValues, const mat covariates) {
    vec placeholder(covariates.rows()) ;
    inputdata dataObject(placeholder, spCoords, timeValues, covariates) ;
    m_predictData = dataObject ;
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
  double m_nuggetSD{ 1e-6 };
  std::string m_distMethod ;

  int m_M{ 0 } ;
  int m_Mlon{ 0 } ;
  int m_Mlat{ 0 } ;
  int m_Mtime{ 0 } ;
  int m_numTips{ 0 } ;
  double m_errorSD{ 0 } ;
  double m_fixedEffSD{ 0 } ;
  bool m_normalHyperprior{ false } ;

  std::unique_ptr<TwoParsProbDist> m_MaternParsHyperparsRhoSpace ;
  std::unique_ptr<TwoParsProbDist> m_MaternParsHyperparsSmoothnessSpace ;
  std::unique_ptr<TwoParsProbDist> m_MaternParsHyperparsRhoTime ;
  std::unique_ptr<TwoParsProbDist> m_MaternParsHyperparsSmoothnessTime ;
  std::unique_ptr<TwoParsProbDist> m_MaternParsHyperparsScaling ;
  std::unique_ptr<TwoParsProbDist> m_ParsHyperparsFixedEffsSD ;
  std::unique_ptr<TwoParsProbDist> m_ParsHyperparsErrorSD ;

  maternVec m_MaternParsSpace ;
  maternVec m_MaternParsTime ;
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
