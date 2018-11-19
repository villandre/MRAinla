#include "TreeNode.h"

namespace MRAinla
{
class TipNode:public TreeNode
{
public:

  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChild(TreeNode *) {assert(false) ;}

  std::vector<TreeNode *> GetChildren() {
    std::vector<TreeNode *> myVec;
    myVec.push_back(NULL) ;
    return myVec;
  } // An input node returns a null pointer when it is asked to provide the address of a child.

  void RemoveChildren() {}
  uint GetM() {return m_depth ;}

  void DeriveAtilde() {
    for (uint k = 0; k < m_depth; k++) {
      for (uint l = 0; l <= k; l++) {
        m_AtildeList.at(k).at(l) = trans(GetB(k)) * m_SigmaInverse * GetB(l) ;
      }
    }
  }

  void DeriveOmega(const inputdata & dataset) {
    std::cout << "Entering DeriveOmega in Tips \n";
    arma::vec subResponses = dataset.responseValues.elem(m_obsInNode) ;

    // std::transform(GetBlist().begin(), std::prev(GetBlist().end()), m_omegaTilde.begin(),
    // [this, subResponses] (arma::mat & Bmatrix) {
    //   std::cout << "Using B... " ;
    //   arma::vec noTemp = arma::trans(Bmatrix) * m_SigmaInverse * subResponses ;
    //   std::cout << "Done! \n" ;
    //   return noTemp;}) ;
    for (uint i = 0 ; i < m_depth; i++) {
      m_omegaTilde.at(i) = arma::trans(GetB(i)) * m_SigmaInverse * subResponses ;
    }
    std::cout << "Leaving DeriveOmega in Tips \n" ;
  }

  void DeriveU(const inputdata & dataset) {
    arma::mat uInMat = arma::trans(dataset.responseValues.elem(m_obsInNode)) * m_SigmaInverse *
      dataset.responseValues.elem(m_obsInNode) ;
    m_u = uInMat.at(0,0) ; // uInMat is supposed to be a 1x1 matrix.
  }

  void DeriveD() {
    std::cout << "Entering DeriveD in tip \n" ;
    double val = 0;
    double sign = 0;
    arma::log_det(val, sign, GetSigma()) ;
    m_d = gsl_sf_log(sign) + val ;
  }

  void ComputeWmat() {
    baseComputeWmat() ;
    m_SigmaInverse = arma::inv_sympd(GetSigma()) ;
  }

  void genRandomKnots(inputdata & dataset, uint & numKnots, const gsl_rng * RNG) {
    m_knotsCoor = spatialcoor(dataset.spatialCoords.rows(m_obsInNode),
                              dataset.timeCoords.elem(m_obsInNode)) ;
  }

  TipNode(const dimensions & dims, const uint & depth, TreeNode * parent,
          const inputdata & dataset, const arma::vec & covPars) {
    baseInitialise(dims, depth, parent, dataset, covPars) ;
    deriveObsInNode(dataset) ;
  }

protected:

  arma::mat GetB(uint & l) {
    if (l == m_depth) {
      throw Rcpp::exception("Trying to get B^M(j_1, ..., j_M)... \n") ;
    }
    return m_Wlist.at(l) ;
  }

  std::vector<arma::mat> GetBlist() {
    return m_Wlist ;
  };

  arma::mat GetSigma() {
    return m_Wlist.at(m_depth) ;
  }

  std::vector<std::vector<arma::mat>> * m_Umat ;
  std::vector<std::vector<arma::mat>> * m_Vmat ;
  std::vector<arma::mat> m_bPred ;
  std::vector<arma::mat> m_bKnots ;
  std::vector<std::vector<arma::mat>> m_Btilde ;
  arma::mat m_SigmaInverse ;
};
}
