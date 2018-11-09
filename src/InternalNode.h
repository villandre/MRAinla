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

  bool IsSolved() {return m_isSolved ;};
  bool CanSolve() ;
  void SetSolved(bool status)  {m_isSolved = status ;};
  InternalNode(dimensions & dims, uint & depth, TreeNode * parent, inputdata & dataset) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    deriveObsInNode(dataset) ;
  }

  InternalNode(dimensions & dims, uint & depth) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = this ;
  }

protected:

   bool m_isSolved ;
   std::vector<TreeNode *> m_children ;
   void genRandomKnots(inputdata &, uint &, const gsl_rng *) ;
};
}
#endif /* INTERMEDIATENODE_H */
