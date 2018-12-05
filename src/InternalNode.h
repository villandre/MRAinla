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
  void DeriveOmega(const inputdata &, const arma::vec &) ;
  void DeriveU(const inputdata &) ;
  void DeriveD() ;
  void ComputeWmat(const arma::vec &) ;
  void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const arma::vec &, const arma::vec &) ;

  void genRandomKnots(inputdata &, uint &, const gsl_rng *) ;

  InternalNode(const dimensions & dims, const uint & depth, TreeNode * parent,
               const inputdata & dataset) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_obsInNode = deriveObsInNode(dataset) ;
    m_omega.resize(m_depth + 1) ;
    m_Alist.resize(m_depth+1) ;
    for (uint i = 0; i < m_Alist.size(); i++) {
      m_Alist.at(i).resize(i+1) ;
    }
  }

  InternalNode(const dimensions & dims, const inputdata & dataset) {
    baseInitialise(dims, 0, this, dataset) ;
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

  // Prediction elements

  void computeEtaTildeParas() ;
  GaussDistParas m_etaTilde ;
};
}
#endif /* INTERMEDIATENODE_H */
