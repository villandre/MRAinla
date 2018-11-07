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
  InternalNode(const vec &, const vec &, const uvec &, const mat &, const uvec &, const uint &) ;

  InternalNode(): _isSolved(false) {
    // TO_DO
  };

protected:

   bool _isSolved ;
   std::vector<TreeNode *> _children ;
};

#endif /* INTERMEDIATENODE_H */
