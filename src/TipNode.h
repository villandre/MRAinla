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
  void SetAlist(const arma::mat &, const uint & uint1, const uint &) { assert(false) ;}
  arma::mat GetAlist(const uint & i, const uint & j) {
    throw Rcpp::exception("Trying to get A matrix at depth M.\n") ;
    arma::mat aMat;
    return aMat;
  }
  void SetKtilde(arma::mat & matrix) { assert(false) ;};
  void SetKtildeInverse(arma::mat & matrix) {assert(false) ;};
  arma::mat GetKtilde() {
    throw Rcpp::exception("Trying to get Ktilde from tip node (depth M).") ;
    return arma::mat(1,1,arma::fill::zeros) ; // This is for the compiler: the throw should prevent the program from getting there.
  }
  arma::mat GetKtildeInverse() {
    throw Rcpp::exception("Trying to get KtildeInverse from tip node (depth M).") ;
    return arma::mat(1,1,arma::fill::zeros) ;
  }
  void DeriveAtilde() {
    for (uint k = 0; k < m_depth+1; k++) {
      for (uint l = 0; l <= k; l++) {
        m_AtildeList.at(k).at(l) = trans(m_Blist.at(k)) * m_SigmaInverse * m_Blist.at(l) ;
      }
    }
  }

  void genRandomKnots(inputdata & dataset, uint & numKnots, const gsl_rng * RNG) {
    m_knotsCoor = spatialcoor(dataset.spatialCoords.rows(m_obsInNode),
                              dataset.timeCoords.elem(m_obsInNode)) ;
  }

  TipNode(dimensions & dims, uint & depth, TreeNode * parent, inputdata & dataset, double & covarianceParameter) {
    baseInitialise(dims, depth, parent, dataset, covarianceParameter) ;
    deriveObsInNode(dataset) ;
  }

protected:
  std::vector<std::vector<arma::mat>> * m_Umat ;
  std::vector<std::vector<arma::mat>> * m_Vmat ;
  std::vector<arma::mat> m_bPred ;
  std::vector<arma::mat> m_bKnots ;
  std::vector<std::vector<arma::mat>> m_Btilde ;
};
}
