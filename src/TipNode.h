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
  void ComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const double & nuggetSD, const std::string & distMethod) {
    baseComputeWmat(covParasSp, covParasTime, scaling, nuggetSD, distMethod) ;

    // m_Wlist.at(m_depth).triangularView<Upper>() = m_Wlist.at(m_depth).triangularView<Lower>() ; // Will this cause aliasing?
    SetKmatrixInverse() ;
    m_K = GetKmatrixInverse().selfadjointView<Eigen::Upper>().ldlt().solve(mat::Identity(GetKmatrixInverse().rows(), GetKmatrixInverse().cols())) ; // The K matrix is some sort of covariance matrix, so it should always be symmetrical..
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

  double GetBelement(const uint & l, const uint & row, const uint & col) {
    if (l == m_depth) {
      return m_Wlist.at(l)(row, m_keepKnotIndices(col)) ;
    } else {
      return m_Wlist.at(l)(row, col) ;
    }
  }

  void SetKmatrixInverse() {
    if (m_knotsThinningRate < 1) {
      mat newMatToAdd = mat(m_keepKnotIndices.size(), m_keepKnotIndices.size()) ;

      for (int i = 0; i < newMatToAdd.rows(); i++) {
        for (int j = 0; j < newMatToAdd.cols(); j++) {
          newMatToAdd(i,j) = m_Wlist.at(m_depth)(m_keepKnotIndices(i), m_keepKnotIndices(j)) ;
        }
      }
      m_Kinverse = newMatToAdd ;
    } else {
      m_Kinverse = m_Wlist.at(m_depth) ;
    }
  };

  mat GetKmatrixInverse() { return m_Kinverse ;}

  // void SetUncorrSD(const double & sd) {
  //   m_uncorrSD = sd ;
  // }
  // mat & GetUpred(const uint & l) { return m_UmatList.at(l) ;}
  double & GetUpredElement(const uint & l, const uint & row, const uint & col) {
    if (l == m_depth) {
      return m_UmatList.at(l)(row, m_keepKnotIndices(col)) ;
    } else {
      return m_UmatList.at(l)(row, col) ;
    }
  }
  std::vector<mat> & GetUmatList() { return m_UmatList ;}

  void SetPredictLocations(const inputdata &) ;
  Eigen::ArrayXi & GetPredIndices() { return m_predsInNode ;}
  int GetNumPreds() override { return m_predsInNode.size() ;}
  void computeUpred(const maternVec &, const maternVec &, const double &, const spatialcoor &, const double &, const std::string &) ;

  void genKnotsOnCube(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & RNG, Eigen::Array<bool, Eigen::Dynamic, 1> &) {

    m_knotsCoor = spatialcoor(rows(dataCoor.spatialCoords, m_obsInNode),
                              elem(dataCoor.timeCoords, m_obsInNode)) ;

    if (m_knotsThinningRate < 1) {
      Eigen::ArrayXi indicesToSampleFrom = Eigen::VectorXi::LinSpaced(
        m_obsInNode.size(),
        0,
        m_obsInNode.size() - 1
      ).array() ;
      m_keepKnotIndices = sampleWithoutReplacement(
        indicesToSampleFrom,
        ceil(double(m_obsInNode.size()) * m_knotsThinningRate),
        true,
        RNG
      ) ;
    }
  }

  void genRandomKnots(spatialcoor & dataCoor, int & numKnots, std::mt19937_64 & RNG) {
    m_knotsCoor = spatialcoor(rows(dataCoor.spatialCoords, m_obsInNode),
                              elem(dataCoor.timeCoords, m_obsInNode)) ;
    if (m_knotsThinningRate < 1) {
      Eigen::ArrayXi indicesToSampleFrom = Eigen::VectorXi::LinSpaced(
        m_obsInNode.size(),
        0,
        m_obsInNode.size() - 1
      ).array() ;
      m_keepKnotIndices = sampleWithoutReplacement(
        indicesToSampleFrom,
        ceil(double(m_obsInNode.size()) * m_knotsThinningRate),
        true,
        RNG
      ) ;
    }
  }

  TipNode(const dimensions & dims, const uint & depth, TreeNode * parent,
          const inputdata & dataset, const double & thinningRate) : m_knotsThinningRate(thinningRate) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_UmatList.resize(m_depth + 1) ;
    m_obsInNode = deriveObsInNode(dataset) ;
    m_keepKnotIndices = Eigen::VectorXi::LinSpaced(GetNumObs(), 0, GetNumObs() - 1).array() ;
  }

  Eigen::ArrayXi GetKeepKnotIndices() override {
    return m_keepKnotIndices ;
  }

protected:

  // double m_uncorrSD{ -1 } ;
  mat m_Kinverse ;

  // Prediction components (should probably be freed once computations are done)

  Eigen::ArrayXi m_predsInNode ;
  std::vector<mat> m_UmatList ;
  Eigen::ArrayXi m_keepKnotIndices ;
  double m_knotsThinningRate{ 1 } ;
};
}
