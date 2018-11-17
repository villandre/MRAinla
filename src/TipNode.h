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
    for (uint k = 0; k < m_depth+1; k++) {
      for (uint l = 0; l <= k; l++) {
        m_AtildeList.at(k).at(l) = trans(m_Blist.at(k)) * m_SigmaInverse * m_Blist.at(l) ;
      }
    }
  }

  void DeriveOmega(const inputdata & dataset) {
    arma::vec subResponses = dataset.responseValues.elem(m_obsInNode) ;
    std::transform(m_Blist.begin(), m_Blist.end(), m_omegaTilde.begin(),
                   [this, subResponses] (arma::mat & Bmatrix) {
                     arma::vec noTemp = arma::trans(Bmatrix) * m_SigmaInverse * subResponses ;
                     return noTemp;}) ;
  }

  void DeriveU(const inputdata & dataset) {
    arma::mat uInMat = arma::trans(dataset.responseValues.elem(m_obsInNode)) * m_SigmaInverse *
      dataset.responseValues.elem(m_obsInNode) ;
    m_u = uInMat.at(0,0) ; // uInMat is supposed to be a 1x1 matrix.
  }

  void DeriveD() {
    double val = 0;
    double sign = 0;
    arma::log_det(val, sign, m_Sigma) ;
    m_d = gsl_sf_log(sign) + val ;
  }
  void DeriveB() {
    std::copy(m_Wlist.begin(), m_Wlist.end(), m_Blist.begin()) ;
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
  std::vector<std::vector<arma::mat>> * m_Umat ;
  std::vector<std::vector<arma::mat>> * m_Vmat ;
  std::vector<arma::mat> m_bPred ;
  std::vector<arma::mat> m_bKnots ;
  std::vector<std::vector<arma::mat>> m_Btilde ;
};
}
