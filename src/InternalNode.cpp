#include "InternalNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

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

  double minTime = (double) min(m_dimensions.time) ;
  uint timeRangeSize = range(m_dimensions.time) ;
  uvec time(numKnots) ;

  time.imbue( [&]() { return gsl_rng_uniform_int(RNG, timeRangeSize) + minTime; } ) ;

  m_knotsCoor = spatialcoor(knotsSp, time) ;
}
