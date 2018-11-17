#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H
namespace MRAinla
{
class InternalNode : public TreeNode
{
public:
  void AddChild(TreeNode * child)  {m_children.push_back(child) ;};
  void RemoveChild(TreeNode *) ;
  std::vector<TreeNode *> GetChildren() {return m_children;};
  void RemoveChildren()  {m_children.clear() ;} ;
  uint GetM() {
    TreeNode * currentAddress = this ;
    while (currentAddress->GetChildren().at(0) != NULL) {
      currentAddress = currentAddress->GetChildren().at(0) ;
    }
    return currentAddress->GetDepth() ;
  }
  void DeriveAtilde() ;
  void DeriveOmega(const inputdata &) ;
  void DeriveU(const inputdata &) ;
  void DeriveD() ;
  void DeriveB() {assert(false) ;}

  void genRandomKnots(inputdata &, uint &, const gsl_rng *) ;

  InternalNode(const dimensions & dims, const uint & depth, TreeNode * parent,
               const inputdata & dataset, const arma::vec & covPars) {
    baseInitialise(dims, depth, parent, dataset, covPars) ;
    deriveObsInNode(dataset) ;
    m_Alist.resize(m_depth+1) ;
    for (uint i = 0; i < m_Alist.size(); i++) {
      m_Alist.at(i).resize(i+1) ;
    }
  }

  InternalNode(const dimensions & dims, const inputdata & dataset, const arma::vec & covPars) {
    baseInitialise(dims, 0, this, dataset, covPars) ;
    uint numObs = dataset.responseValues.size() ;
    m_obsInNode = arma::regspace<arma::uvec>(0, numObs - 1) ;
    m_omega.resize(m_depth + 1) ;
    m_Alist.resize(m_depth+1) ;
    for (uint i = 0; i < m_Alist.size(); i++) {
      m_Alist.at(i).resize(i+1) ;
    }
  }

protected:

  std::vector<TreeNode *> m_children ;
  std::vector<std::vector<arma::mat>> m_Alist ;
  std::vector<arma::vec> m_omega ;
  arma::mat m_Ktilde ;
  arma::mat m_KtildeInverse ;
};
}
#endif /* INTERMEDIATENODE_H */
