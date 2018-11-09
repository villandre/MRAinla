#include "InternalNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

bool InternalNode::CanSolve()
{
  std::vector<bool> childDefined(m_children.size()) ;
  childDefined.reserve(m_children.size()) ;
  for (auto & i : m_children)
  {
    childDefined.push_back(i->IsSolved()) ;
  }
  return std::all_of(childDefined.begin(), childDefined.end(), [](bool v) { return v; });
}

void InternalNode::RemoveChild(TreeNode * childToRemove)
{
  auto ChildIterPos = std::find(m_children.begin(), m_children.end(), childToRemove) ;
  if (ChildIterPos == m_children.end())
  {
    cerr << "Warning: Trying to remove a child that was not found! \n" ;
  }
  else
  {
    m_children.erase(ChildIterPos) ;
  }
}

void InternalNode::genRandomKnots(inputdata & dataset, uint & numKnots, const gsl_rng * RNG) {

  mat knotsSp(numKnots, 2) ;
  double minLon = (double) min(m_dimensions.longitude) ;
  double maxLon = (double) max(m_dimensions.longitude) ;
  double minLat = (double) min(m_dimensions.latitude) ;
  double maxLat = (double) max(m_dimensions.latitude) ;

  for (mat::iterator iter = knotsSp.begin() ; iter != std::prev(knotsSp.end()) ; std::advance(iter, 2)) {
    (*iter) = gsl_ran_flat(RNG, minLon, maxLon) ;
    (*std::next(iter)) = gsl_ran_flat(RNG, minLat, maxLat) ;
  }
}
