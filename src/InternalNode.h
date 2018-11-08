#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H

class InternalNode : public TreeNode
{
public:
  void AddChild(TreeNode * child) {_children.push_back(child) ;};
  void RemoveChild(TreeNode *) ;
  std::vector<TreeNode *> GetChildren() {return _children;};
  void RemoveChildren() {_children.clear() ;} ;

  bool IsSolved() {return _isSolved ;};
  bool CanSolve() ;
  void SetSolved(bool status) {_isSolved = status ;};
  InternalNode(dimtype & dims, uint & depth, TreeNode * parent, datasettype & dataset) {
    _dimensions = dims;
    _depth = depth ;
    _parent = parent ;
    deriveObsInNode(dataset) ;
  }

  InternalNode(dimtype & dims, uint & depth) {
    _dimensions = dims;
    _depth = depth ;
    _parent = this ;
  }

protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
   void genRandomKnots(datasettype & dataset) {
     // _knotsCoor = std::make_tuple( , ) ; //TO_DO
   }
};

#endif /* INTERMEDIATENODE_H */
