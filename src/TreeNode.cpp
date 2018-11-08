#include "TreeNode.h"

void TreeNode::deriveObsInNode(datasettype & dataset) {
  uvec lonCheck = (std::get<1>(dataset).col(0) > std::get<0>(_dimensions).at(0)) *
    (std::get<1>(dataset).col(0) <= std::get<0>(_dimensions).at(1)) ; // Longitude check
  uvec latCheck = (std::get<1>(dataset).col(1) > std::get<1>(_dimensions).at(0)) *
    (std::get<1>(dataset).col(1) <= std::get<1>(_dimensions).at(1)) ; // Latitude check
  uvec timeCheck = (std::get<2>(dataset) > std::get<2>(_dimensions).at(0)) *
    (std::get<2>(dataset) <= std::get<2>(_dimensions).at(1)) ; // Time check
  _obsInNode = find(lonCheck * latCheck * timeCheck) ;
}
