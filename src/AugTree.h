#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

namespace MRAinla {

class AugTree
{
public:
  AugTree(uint &, arma::vec &, arma::vec &, arma::uvec &, arma::vec &, arma::mat &, arma::uvec &, uint &, unsigned long int &, arma::mat &) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

  void ComputeMRAloglik() ;
  double ComputeGlobalLogLik() ;
  void ComputeConditionalPrediction(const spatialcoor &) ;
  arma::mat ComputePosteriors(spatialcoor &, double &) ;

  double GetLoglik() {return m_logLik ;}
  gsl_rng * GetRandomNumGenerator() {return m_randomNumGenerator ;}
  inputdata GetDataset() {return m_dataset;}
  uint GetNumTips() {return m_numTips ;}

  void SetRNG(gsl_rng * myRNG) { m_randomNumGenerator = myRNG ;}

  void SetCovParameters(arma::vec & covParas) {m_covParameters = covParas ;}
  void SetFixedEffParameters(arma::vec & fixedParas) {m_fixedEffParameters = fixedParas ;}
  void SetFixedEffSD(const double & fixedEffSD) {m_fixedEffSD = fixedEffSD ;}
  void SetErrorSD(const double & errorSD) {m_errorSD = errorSD ;}

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

  arma::vec m_covParameters ;
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
  void computeOmegas() ;

  void computeU() ;
  void computeD() ;

  // Prediction functions

  std::vector<TreeNode *> GetLevel(const uint) ;
  void distributePredictionData(const spatialcoor &) ;
  void computeBtildeInTips() ;
  spatialcoor m_predictLocations ;
};
}
#endif
