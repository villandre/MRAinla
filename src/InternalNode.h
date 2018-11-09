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

  void genRandomKnots(inputdata &, uint &, const gsl_rng *) ;

  InternalNode(dimensions & dims, uint & depth, TreeNode * parent, inputdata & dataset) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    deriveObsInNode(dataset) ;
  }

  InternalNode(dimensions & dims, inputdata & dataset) {
    m_dimensions = dims;
    m_depth = 0 ;
    uint numObs = dataset.responseValues.size() ;
    m_obsInNode = arma::regspace<arma::uvec>(0, numObs - 1) ;
    m_parent = this ;
  }

protected:

   std::vector<TreeNode *> m_children ;
};
}
#endif /* INTERMEDIATENODE_H */
