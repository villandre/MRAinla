#include "TreeNode.h"

class TipNode:public TreeNode
{
public:

  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChild(TreeNode *) {assert(false) ;}
  std::vector<TreeNode *> GetChildren() {std::vector<TreeNode *> myVec; myVec.push_back(NULL) ; return myVec;} // An input node returns a null pointer when it is asked to provide the address of a child.
  void RemoveChildren() {}

  bool IsSolved() {return true ;}
  bool CanSolve() {return true ;}
  void SetSolved(bool status) {}

  TipNode()
  {
    // TO_DO
  }

protected:
  std::vector<std::vector<mat>> * _Umat ;
  std::vector<std::vector<mat>> * _Vmat ;
  std::vector<mat> _bPred ;
  std::vector<mat> _bKnots ;
  std::vector<std::vector<mat>> _Btilde ;
};
