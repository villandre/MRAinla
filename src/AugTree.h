#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

namespace MRAinla {

class AugTree
{
public:
  AugTree(uint &, arma::vec &, arma::vec &, arma::uvec &, arma::vec &, arma::mat &, arma::uvec &, uint &, unsigned long int &, arma::mat &) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

  double ComputeMRAlogLik(const arma::vec &, const bool) ;
  std::vector<GaussDistParas> ComputeConditionalPrediction(const inputdata &) ;
  arma::mat ComputePosteriors(spatialcoor &, double &) ;
  double ComputeLogJointCondTheta(const arma::vec &, const arma::vec &, const arma::vec &, const double, const double) ;
  double ComputeGlobalLogLik(const arma::vec &, const arma::vec &, const double) ;
  double ComputeLogFullConditional(const arma::vec &, const arma::vec &) ;
  double ComputeLogPriors(const arma::vec &, const double, const double, const double, const double) ;

  // double GetMRAlogLik() const {return m_MRAlogLik ;}
  gsl_rng * GetRandomNumGenerator() {return m_randomNumGenerator ;}
  inputdata GetDataset() {return m_dataset;}
  uint GetNumTips() {return m_numTips ;}
  // arma::vec GetCovParameters() {return m_covParameters ;}

  void SetRNG(gsl_rng * myRNG) { m_randomNumGenerator = myRNG ;}

  // void SetCovParameters(arma::vec & covParas) {m_covParameters = covParas ;}
  void SetFixedEffParameters(arma::vec & fixedParas) {m_fixedEffParameters = fixedParas ;}
  void SetFixedEffSD(const double & fixedEffSD) {m_fixedEffSD = fixedEffSD ;}
  void SetErrorSD(const double & errorSD) {m_errorSD = errorSD ;}
  void CleanPredictionComponents() ;
  void CenterResponse() ;
  arma::sp_mat createHstar() ;
  arma::sp_mat createSigmaStarInverse() ;
  arma::sp_mat createQ() ;

  double ComputeJointPsiMarginal(const arma::vec &, const double, const double, const double, const double) ;

  ~ AugTree() {
    deallocate_container(m_vertexVector) ;
    gsl_rng_free(m_randomNumGenerator) ;};

private:

  std::vector<TreeNode *> m_vertexVector ;
  std::vector<TreeNode *> getLevelNodes(uint & level) ;

  double m_MRAlogLik{ 0 } ;
  double m_globalLogLik{ } ;

  uint m_M{ 0 } ;
  uint m_numTips{ 0 } ;
  double m_errorSD ;
  double m_fixedEffSD ;
  arma::vec m_spatialComponents ;

  arma::vec m_MRAcovParas ;
  arma::vec m_fixedEffParameters ;

  dimensions m_mapDimensions;

  inputdata m_dataset ; // First element is response, second is spatial coordinates, last is time.
  // Note that the generator's seed is determined by the system variable
  // GSL_RNG_SEED and takes value 0 by default.
  gsl_rng * m_randomNumGenerator ;

  // Tree construction functions //
  void BuildTree(const uint &) ;
  void createLevels(TreeNode *, const uint &) ;
  void generateKnots(TreeNode *) ;

  // Likelihood computations functions

  void computeWmats() ;
  void deriveAtildeMatrices() ;
  void computeOmegas(const arma::vec &) ;

  void computeU(const arma::vec &) ;
  void computeD() ;

  // Prediction functions

  std::vector<TreeNode *> GetLevel(const uint) ;
  void distributePredictionData(const spatialcoor &) ;
  void computeBtildeInTips() ;
  spatialcoor m_predictLocations ;

  // INLA functions
  std::vector<arma::mat *> getKmatricesInversePointers() {
    std::vector<arma::mat *> KmatrixInverseList ;
    for (auto & i : m_vertexVector) KmatrixInverseList.push_back(i->GetKmatrixInverseAddress()) ;
    return KmatrixInverseList ;
  }
  arma::sp_mat createFmatrix() ;
};
}
#endif
