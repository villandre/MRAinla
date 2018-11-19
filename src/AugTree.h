#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

namespace MRAinla {

class AugTree
{
public:
  AugTree(uint &, arma::vec &, arma::vec &, arma::uvec &, arma::vec &, arma::mat &, arma::uvec &, uint &, unsigned long int &, arma::vec & covPars) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

  void ComputeLoglik() ;

  double GetLoglik() {return m_logLik ;}
  gsl_rng * GetRandomNumGenerator() {return m_randomNumGenerator ;}
  inputdata GetDataset() {return m_dataset;}
  uint GetNumTips() {return m_numTips ;}

  void SetRNG(gsl_rng * myRNG) { m_randomNumGenerator = myRNG ;}

  void ComputeLoglik(Rcpp::List &, Rcpp::List &, Rcpp::NumericVector &) ;

  ~ AugTree() {
    deallocate_container(m_vertexVector) ;
    gsl_rng_free(m_randomNumGenerator) ;};

private:

  std::vector<TreeNode *> m_vertexVector ;
  std::vector<TreeNode *> getLevelNodes(uint & level) ;

  double m_logLik{ 0 } ;

  uint m_M{ 0 } ;
  uint m_numTips{ 0 } ;

  dimensions m_mapDimensions;

  inputdata m_dataset ; // First element is response, second is spatial coordinates, last is time.
  // Note that the generator's seed is determined by the system variable
  // GSL_RNG_SEED and takes value 0 by default.
  gsl_rng * m_randomNumGenerator ;

  // Tree construction functions //
  void BuildTree(const uint &, const arma::vec &) ;
  void createLevels(TreeNode *, const uint &, const arma::vec &) ;
  void generateKnots(TreeNode *) ;

  // Likelihood computations functions

  void computeWmats() ;
  void deriveAtildeMatrices() ;
  void computeOmegas() ;

  void computeU() ;
  void computeD() ;
};
}
#endif
