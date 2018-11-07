#include "InternalNode.h"

void InternalNode::InvalidateSolution()
{
  _isSolved = false ;

  if (_parent != NULL)
  { // Root has a NULL parent.
    _parent->InvalidateSolution() ;
  }
}

bool InternalNode::CanSolve()
{
  std::vector<bool> childDefined(_children.size()) ;
  childDefined.reserve(_children.size()) ;
  for (auto & i : _children)
  {
    childDefined.push_back(i->IsSolved()) ;
  }
  return std::all_of(childDefined.begin(), childDefined.end(), [](bool v) { return v; });
}

void InternalNode::ComputeSolutions()
{
  // TO_DO
}

void InternalNode::ComputeSolution()
{
  // TO_DO
}

void InternalNode::RemoveChild(TreeNode * childToRemove)
{
  auto ChildIterPos = std::find(_children.begin(), _children.end(), childToRemove) ;
  if (ChildIterPos == _children.end())
  {
    cerr << "Warning: Trying to remove a child that was not found! \n" ;
  }
  else
  {
    _children.erase(ChildIterPos) ;
  }
}
