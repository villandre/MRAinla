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
  cout << "Number of knots: " << numKnots << "\n \n" ;

  double minLon = min(m_dimensions.longitude) ;
  double maxLon = max(m_dimensions.longitude) ;

  double minLat = min(m_dimensions.latitude) ;
  double maxLat = max(m_dimensions.latitude) ;
  cout << "Generating coordinates for lon/lat... " ;
  for (mat::iterator iter = knotsSp.begin() ; iter != std::prev(knotsSp.end()) ; std::advance(iter, 2)) {
    (*iter) = gsl_ran_flat(RNG, minLon, maxLon) ;
    (*std::next(iter)) = gsl_ran_flat(RNG, minLat, maxLat) ;
  }
  cout << "Done! \n \n" ;

  double minTime = (double) min(m_dimensions.time) ;
  uint timeRangeSize = range(m_dimensions.time) ;
  uvec time(numKnots) ;
  cout << "Creating time... " ;
  time.imbue( [&]() { return gsl_rng_uniform_int(RNG, timeRangeSize) + minTime; } ) ;
  cout << "Done! \n \n" ;
  m_knotsCoor = spatialcoor(knotsSp, time) ;
}
