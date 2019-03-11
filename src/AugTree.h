#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

struct IGhyperParas{
  // Default IG mean is 1 (mean = beta/(alpha - 1)), var is beta^2/[(alpha-1)^2(alpha-2)]
  // If alpha <= 2, variance is infinite,
  double m_alpha = 3;
  double m_beta = 2;
  IGhyperParas() {}
  IGhyperParas(double alpha, double beta) : m_alpha(alpha), m_beta(beta) {
    // if (alpha < 2) {
    //   std::cout << "Careful: IG distribution for hyperparameter has infinite variance! \n" ;
    // }
  }
};

namespace MRAinla {

class AugTree
{
public:
  AugTree(uint &, arma::fvec &, arma::fvec &, arma::fvec &, arma::vec &, arma::fmat &, arma::fvec &, uint &, unsigned long int &, arma::fmat &) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

  void ComputeMRAlogLik(const bool WmatsAvailable = false) ;
  void ComputeMRAlogLikAlt(const bool WmatsAvailable = false) ;
  std::vector<GaussDistParas> ComputeConditionalPrediction(const inputdata &) ;
  void ComputeLogJointCondTheta() ;
  void ComputeGlobalLogLik() ;
  void ComputeLogFullConditional() ;
  void ComputeLogPriors() ;

  // double GetMRAlogLik() const {return m_MRAlogLik ;}
  gsl_rng * GetRandomNumGenerator() {return m_randomNumGenerator ;}
  inputdata GetDataset() {return m_dataset;}
  uint GetNumTips() {return m_numTips ;}
  // arma::vec GetCovParameters() {return m_covParameters ;}
  std::vector<IGhyperParas> GetMRAcovParasIGalphaBeta() { return m_MRAcovParasIGalphaBeta ;}

  void SetRNG(gsl_rng * myRNG) { m_randomNumGenerator = myRNG ;}

  // void SetCovParameters(arma::vec & covParas) {m_covParameters = covParas ;}
  void SetFixedEffParameters(arma::vec & fixedParas) {
    if (arma::approx_equal(fixedParas, m_fixedEffParameters, "absdiff", 1e-8)) m_recomputeGlobalLogLik = true ;
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
  void SetMRAcovParas(const arma::vec & MRAcovParas) {
    if (arma::approx_equal(MRAcovParas, m_MRAcovParas, "absdiff", 1e-8)) m_recomputeMRAlogLik = true ;
    m_MRAcovParas = MRAcovParas ;
  }
  void SetMRAcovParasIGalphaBeta(const std::vector<IGhyperParas> & alphaBeta) {m_MRAcovParasIGalphaBeta = alphaBeta ;}
  void SetErrorIGalphaBeta(const IGhyperParas & alphaBeta) {m_errorIGalphaBeta = alphaBeta ;}
  void SetFixedEffIGalphaBeta(const IGhyperParas & alphaBeta) {m_fixedEffIGalphaBeta = alphaBeta ;}
  void SetFEmu(const arma::vec & muVec) {m_FEmu = muVec ;}

  void CleanPredictionComponents() ;
  void CenterResponse() ;
  arma::sp_mat createHstar() ;
  arma::sp_mat CombineKandFEmatrices() ;
  arma::sp_mat createQ() ;

  double ComputeLogJointPsiMarginal() ;
  // double ComputeJointPsiMarginalPropConstant(const arma::vec &, const double, const double, const double, const double) ;

  ~ AugTree() {
    deallocate_container(m_vertexVector) ;
    gsl_rng_free(m_randomNumGenerator) ;};

private:

  std::vector<TreeNode *> m_vertexVector ;
  std::vector<TreeNode *> getLevelNodes(uint & level) ;
  bool m_recomputeMRAlogLik{ true } ; // When this flag is true, more computations are required to get the log-lik.
  bool m_recomputeGlobalLogLik{ true } ; // When this flag is true, the global log-likelihood (conditional on all mean parameters) needs to be recomputed.

  double m_MRAlogLik{ 0 } ;
  double m_globalLogLik{ 0 } ;
  double m_fixedEffContrib{ 0 } ;
  double m_logPrior{ 0 } ;
  double m_logCondDist{ 0 } ;
  double m_logFullCond{ 0 } ;
  uint m_numKnots ;

  uint m_M{ 0 } ;
  uint m_numTips{ 0 } ;
  double m_errorSD ;
  double m_fixedEffSD ;
  arma::vec m_spatialComponents ;
  std::vector<IGhyperParas> m_MRAcovParasIGalphaBeta ;
  IGhyperParas m_fixedEffIGalphaBeta ;
  IGhyperParas m_errorIGalphaBeta ;

  arma::vec m_MRAcovParas ;
  arma::vec m_fixedEffParameters ;
  arma::vec m_FEmu ;

  dimensions m_mapDimensions;

  inputdata m_dataset ; // First element is response, second is spatial coordinates, last is time.
  // Note that the generator's seed is determined by the system variable
  // GSL_RNG_SEED and takes value 0 by default.
  gsl_rng * m_randomNumGenerator ;

  std::vector<TreeNode *> Descendants(std::vector<TreeNode *>) ;
  void diveAndUpdate(TreeNode *, std::vector<TreeNode *> *) ;

  // Tree construction functions //
  void BuildTree(const uint &) ;
  void createLevels(TreeNode *, const uint &) ;
  void generateKnots(TreeNode *) ;

  // Likelihood computations functions

  void computeWmats() ;
  void deriveAtildeMatrices() ;
  void computeOmegas() ;

  void computeU() ;
  void computeD() ;

  // Prediction functions

  std::vector<TreeNode *> GetLevel(const uint) ;
  void distributePredictionData(const spatialcoor &) ;
  void computeBtildeInTips() ;
  spatialcoor m_predictLocations ;

  // INLA functions
  arma::vec m_MRArandomValues ;
  arma::vec m_MRAetaValues ;
  arma::vec m_Vstar ;
  arma::vec m_MRAvalues ;
  std::vector<arma::mat *> getKmatricesInversePointers() {
    std::vector<arma::mat *> KmatrixInverseList ;
    for (auto & i : m_vertexVector) KmatrixInverseList.push_back(i->GetKmatrixInverseAddress()) ;
    return KmatrixInverseList ;
  }
  arma::sp_mat createFmatrix() ;
  arma::vec optimJointHyperMarg(const arma::vec &, const double, const double, const double, const double) ;
  double logDeterminantQmat(const arma::sp_mat & Qmat) ;
  arma::uvec extractBlockIndicesFromLowerRight(const arma::sp_mat &) ;
  // void invFromDecomposition(const arma::sp_mat &, const arma::sp_mat &, const arma::sp_mat &, arma::sp_mat *,
  //                                const arma::uvec &) ;
  double logDeterminantFullConditional(const arma::sp_mat &) ;
  arma::vec ComputeFullConditionalMean(const arma::vec &, const arma::sp_mat &) ;
  template <typename MatrixType>
  inline typename MatrixType::Scalar logdet(const MatrixType& M, bool use_cholesky = false) ;

};
}
#endif
