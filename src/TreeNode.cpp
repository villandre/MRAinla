#include "TreeNode.h"

using namespace arma ;
using namespace MRAinla ;

void TreeNode::deriveObsInNode(inputdata & dataset) {
  uvec lonCheck = (dataset.spatialCoords.col(0) > min(m_dimensions.longitude)) *
    (dataset.spatialCoords.col(0) <= max(m_dimensions.longitude)) ; // Longitude check
  uvec latCheck = (dataset.spatialCoords.col(1) > min(m_dimensions.latitude)) *
    (dataset.spatialCoords.col(1) <= max(m_dimensions.latitude)) ; // Latitude check
  uvec timeCheck = (dataset.timeCoords > min(m_dimensions.time)) *
    (dataset.timeCoords <= max(m_dimensions.time)) ; // Time check
  m_obsInNode = find(lonCheck * latCheck * timeCheck) ; // find is equivalent to which in R
}
