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
  void SetAlist(const arma::mat & matrix, const uint & k, const uint & l) {
    m_Alist.at(k).at(l) = matrix ;
  ;}
  arma::mat GetAlist(const uint & i, const uint & j) {return m_Alist.at(i).at(j) ;}
  void SetKtilde(arma::mat & matrix) { m_Ktilde = matrix ;};
  void SetKtildeInverse(arma::mat & matrix) {m_KtildeInverse = matrix ;};
  arma::mat GetKtilde() {return m_Ktilde ;}
  arma::mat GetKtildeInverse() {return m_KtildeInverse ;}
  void DeriveAtilde() ;

  void genRandomKnots(inputdata &, uint &, const gsl_rng *) ;

  InternalNode(dimensions & dims, uint & depth, TreeNode * parent, inputdata & dataset, double & covarianceParameter) {
    baseInitialise(dims, depth, parent, dataset, covarianceParameter) ;
    deriveObsInNode(dataset) ;
  }

  InternalNode(dimensions & dims, inputdata & dataset, double & covarianceParameter) {
    baseInitialise(dims, 0, this, dataset, covarianceParameter) ;
    uint numObs = dataset.responseValues.size() ;
    m_obsInNode = arma::regspace<arma::uvec>(0, numObs - 1) ;
  }

protected:

  std::vector<TreeNode *> m_children ;
  std::vector<std::vector<arma::mat>> m_Alist ;
  arma::mat m_Ktilde ;
  arma::mat m_KtildeInverse ;
};
}
#endif /* INTERMEDIATENODE_H */
