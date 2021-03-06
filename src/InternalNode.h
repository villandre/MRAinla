#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H
namespace MRAinla
{

struct intCube{
  int edgeLengthInUnits ;
  Eigen::Array3d scalingFactors ;
  Eigen::Array3d offsetCoords ;
  intCube() { } ;
  intCube(int edgeLength, Eigen::Array3d scalingFactors, Eigen::Array3d offsetCoords): edgeLengthInUnits(edgeLength),
    scalingFactors(scalingFactors), offsetCoords(offsetCoords) {  }
  std::array<Eigen::Array3d, 8> getCorners(bool units) {
    Eigen::Array3d coef = Eigen::Array3d::Constant(edgeLengthInUnits) ;
    Eigen::Array3d offset = Eigen::Array3d::Zero(edgeLengthInUnits) ;
    if (!units) {
      coef = scalingFactors * (edgeLengthInUnits - 1) ;
      offset = offsetCoords ;
    }
    uint index = 0 ;
    std::array<Eigen::Array3d, 8> corners ;
    for (uint width = 0 ; width < 2 ; width++) {
      for (uint height = 0 ; height < 2 ; height++) {
        for (uint depth = 0 ; depth < 2 ; depth++) {
          corners.at(index) = {width * coef(0) + offset(0),
                               height * coef(1) + offset(1),
                               depth * coef(2) + offset(2) } ;
          index += 1 ;
        }
      }
    }
    return corners ;
  }
  Eigen::Array3d getEdgePointCoor(const uint edgeIndex, const uint pointIndex) {
    Eigen::Array3i pointReturn ;
    switch (edgeIndex) {
    case 0:
      pointReturn  << pointIndex, 0, 0 ;
      break ;
    case 1:
      pointReturn  << (edgeLengthInUnits - 1), 0, pointIndex ;
      break ;
    case 2:
      pointReturn << (edgeLengthInUnits - pointIndex - 1), 0, (edgeLengthInUnits - 1) ;
      break ;
    case 3:
      pointReturn << 0, 0, (edgeLengthInUnits - pointIndex - 1) ;
      break ;
    case 4:
      pointReturn  << pointIndex, (edgeLengthInUnits - 1), 0  ;
      break ;
    case 5:
      pointReturn  << (edgeLengthInUnits - 1), (edgeLengthInUnits - 1), pointIndex ;
      break ;
    case 6:
      pointReturn <<  (edgeLengthInUnits - pointIndex - 1), (edgeLengthInUnits - 1), (edgeLengthInUnits - 1) ;
      break ;
    case 7:
      pointReturn << 0, (edgeLengthInUnits - 1), (edgeLengthInUnits - pointIndex - 1) ;
      break ;
    case 8:
      pointReturn << 0, pointIndex, 0 ;
      break ;
    case 9:
      pointReturn << (edgeLengthInUnits - 1), (edgeLengthInUnits - pointIndex - 1), 0 ;
      break ;
    case 10:
      pointReturn << (edgeLengthInUnits - 1), pointIndex, (edgeLengthInUnits - 1);
      break ;
    case 11:
      pointReturn << 0, (edgeLengthInUnits - pointIndex - 1), (edgeLengthInUnits - 1) ;
      break ;
    } ;
    Eigen::ArrayXd pointReturnRecast = pointReturn.cast<double>() ;
    return pointReturnRecast * scalingFactors + offsetCoords ;
  }
};

class InternalNode : public TreeNode
{
public:
  void AddChild(TreeNode * child)  {m_children.push_back(child) ;};
  void RemoveChild(TreeNode *) ;
  std::vector<TreeNode *> GetChildren() {return m_children;};
  void RemoveChildren()  {m_children.clear() ;} ;
  int GetM() {
    TreeNode * currentAddress = this ;
    while (currentAddress->GetChildren().at(0) != NULL) {
      currentAddress = currentAddress->GetChildren().at(0) ;
    }
    return currentAddress->GetDepth() ;
  }
  void ComputeWmat(const maternVec &, const maternVec &, const double &, const double &, const std::string &) ;

  void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const vec &) ;

  mat & GetB(const uint & l) { throw Rcpp::exception("Trying to get B matrix in internal node.\n") ;}
  double GetBelement(const uint & l, const uint & row, const uint & col) { throw Rcpp::exception("Trying to get elements of B matrix from internal node.\n") ;}
  // void SetUncorrSD(const double &) {throw Rcpp::exception("Trying to add uncorrelated error for internal nodes! \n") ;}
  // mat & GetUpred(const uint & l) { throw Rcpp::exception("Upred matrices only computed for tip nodes! \n") ;}
  double & GetUpredElement(const uint & l, const uint & row, const uint & col) { throw Rcpp::exception("Upred matrices only computed for tip nodes! \n") ; }
  std::vector<mat> & GetUmatList() { throw Rcpp::exception("UmatList only available in tip nodes! \n") ;}
  void SetPredictLocations(const inputdata & data) { throw Rcpp::exception("Trying to attach predict locations to internal nodes! Predict locations should only be defined in the tips! \n") ;}
  Eigen::ArrayXi & GetPredIndices() { throw Rcpp::exception("Prediction locations not defined in internal nodes! \n");}
  void computeUpred(const maternVec &, const maternVec &, const double &, const spatialcoor &, const double &, const std::string &) {
    throw Rcpp::exception("Upred matrices need not be computed in internal nodes! \n") ;
  }
  mat GetKmatrixInverse() {
    return m_Wlist.at(m_depth) ;
  };

  void genKnotsOnCube(spatialcoor &, int &, std::mt19937_64 &, Eigen::Array<bool, Eigen::Dynamic, 1> &) ;
  void genRandomKnots(spatialcoor &, int &, std::mt19937_64 &) ;

  InternalNode(const dimensions & dims, const uint & depth, TreeNode * parent,
               const inputdata & dataset) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_obsInNode = deriveObsInNode(dataset) ;
  }

  InternalNode(const dimensions & dims, const inputdata & dataset) {
    baseInitialise(dims, 0, this, dataset) ;
    int numObs = dataset.responseValues.size() ;
    m_obsInNode = Eigen::VectorXi::LinSpaced(numObs, 0, numObs - 1) ;
  }

protected:

  std::vector<TreeNode *> m_children ;
};
}
#endif /* INTERMEDIATENODE_H */
