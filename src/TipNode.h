#include "TreeNode.h"

namespace MRAinla
{
class TipNode:public TreeNode
{
public:

  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChild(TreeNode *) {assert(false) ;}

  std::vector<TreeNode *> GetChildren() {
    std::vector<TreeNode *> myVec;
    myVec.push_back(NULL) ;
    return myVec;
  } // An input node returns a null pointer when it is asked to provide the address of a child.

  void RemoveChildren() {}
  int GetM() {return m_depth ;}

  void ComputeWmat(const maternVec & covParasSp, const double & scaling, const double & spaceNuggetSD, const std::string & distMethod) {
    baseComputeWmat(covParasSp, scaling, spaceNuggetSD, distMethod) ;
    m_SigmaInverse = GetSigma().selfadjointView<Eigen::Upper>().ldlt().solve(Eigen::MatrixXd::Identity(GetSigma().rows(), GetSigma().cols())) ;
  }

  mat & GetB(const uint & l) {
    if (l > m_depth) {
      Rprintf("Error message... \n") ;
      Rprintf("Error occured in node %i. \n", m_nodeId) ;
      Rprintf("Trying to get index %i while tree has %i layers. \n", l, m_depth) ;
      throw Rcpp::exception("Trying to get B^l(j_1, ..., j_M) with l > M! \n") ;
    }
    return m_Wlist.at(l) ; // This works for l = M because knot positions in tips correspond to observation positions. It wouldn't be valid otherwise.
  }

  mat GetSigma() {
    Eigen::MatrixXd Sigma = m_Wlist.at(m_depth).selfadjointView<Eigen::Upper>() ;
    return Sigma + std::pow(m_uncorrSD, 2) * Eigen::MatrixXd::Identity(m_Wlist.at(m_depth).rows(), m_Wlist.at(m_depth).cols()) ;
  }
  // The following Kmatrix-related functions work because of the correspondence between knots
  // and observations at the finest resolution.
  mat & GetKmatrix() {return m_SigmaInverse ;}
  mat * GetKmatrixAddress() {return &m_SigmaInverse ;}
  mat * GetKmatrixInverseAddress() { return &m_Wlist.at(m_depth) ;}
  mat GetKmatrixInverse() {return GetSigma() - std::pow(m_uncorrSD, 2) * Eigen::MatrixXd::Identity(GetSigma().rows(), GetSigma().cols());}
  void SetUncorrSD(const double & sd) {
    m_uncorrSD = sd ;
  }
  mat & GetUpred(const uint & l) { return m_UmatList.at(l) ;}
  std::vector<mat> & GetUmatList() { return m_UmatList ;}

  void SetPredictLocations(const inputdata &) ;
  Eigen::ArrayXi & GetPredIndices() { return m_predsInNode ;}
  void computeUpred(const maternVec &, const double &, const spatialcoor &, const double &, const std::string &) ;

  void genKnotsOnSquare(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & RNG, Eigen::Array<bool, Eigen::Dynamic, 1> &) {
    m_knotsCoor = spatialcoor(rows(dataCoor.spatialCoords, m_obsInNode)) ;
  }
  void genRandomKnots(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & generator) {
    m_knotsCoor = spatialcoor(rows(dataCoor.spatialCoords, m_obsInNode)) ;
  }

  TipNode(const spaceDimensions & dims, const uint & depth, TreeNode * parent,
          const inputdata & dataset) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_UmatList.resize(m_depth + 1) ;
    m_obsInNode = deriveObsInNode(dataset) ;
  }

protected:

  mat m_SigmaInverse ;
  double m_uncorrSD{ -1 } ;

  // Prediction components (should probably be freed once computations are done)

  Eigen::ArrayXi m_predsInNode ;
  std::vector<mat> m_UmatList ;
};
}
