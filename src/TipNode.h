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

  void DeriveAtilde() {
    for (uint k = 0; k < m_depth; k++) {
      for (uint l = 0; l <= k; l++) {
        m_AtildeList.at(k).at(l) = GetB(k).transpose() * m_SigmaInverse * GetB(l) ;
      }
    }
  }

  void DeriveOmega(const vec &) ;

  void DeriveU(const vec & responseValues) {
    vec subResponse = elem(responseValues.array(), m_obsInNode) ;
    mat uInMat = subResponse.transpose() * m_SigmaInverse *
      subResponse ;
    m_u = uInMat(0,0) ; // uInMat is supposed to be a 1x1 matrix.
  }

  void DeriveD() {
    double val = 0;
    val = GetSigma().colPivHouseholderQr().logAbsDeterminant() ;
    m_d = val ;
  }

  void ComputeWmat(const maternVec & covParasSp, const maternVec & covParasTime, const double & scaling, const bool matern, const double & spaceNuggetSD, const double & timeNuggetSD) {
    baseComputeWmat(covParasSp, covParasTime, scaling, matern, spaceNuggetSD, timeNuggetSD) ;
    // m_Wlist.at(m_depth).triangularView<Eigen::Upper>() = m_Wlist.at(m_depth).triangularView<Eigen::Lower>().eval() ; // Will this cause aliasing?
    m_SigmaInverse = GetSigma().selfadjointView<Eigen::Lower>().ldlt().solve(Eigen::MatrixXd::Identity(GetSigma().rows(), GetSigma().cols())) ;
  }

  // void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const vec &) ;
  std::vector<std::vector<mat>> GetAlist() const {throw Rcpp::exception("Trying to get an A matrix in a tip node! \n") ;}
  mat GetKtilde() const {throw Rcpp::exception("Trying to get Ktilde in a tip node! \n") ;}
  // void deriveBtilde(const spatialcoor & ) ;
  // void computeBpred(const spatialcoor &, const vec &) ;
  // GaussDistParas CombineEtaDelta(const inputdata &, const vec &) ;
  // GaussDistParas GetEtaDelta() const { return m_deltaTilde ;}

  mat GetB(const uint & l) {
    if (l > m_depth) {
      printf("Error message... \n") ;
      printf("Error occured in node %i. \n", m_nodeId) ;
      printf("Trying to get index %i while tree has %i layers. \n", l, m_depth) ;
      throw Rcpp::exception("Trying to get B^l(j_1, ..., j_M) with l > M! \n") ;
    }
    return m_Wlist.at(l) ; // This works for l = M because knot positions in tips correspond to observation positions. It wouldn't be valid otherwise.
  }

  mat GetSigma() {
    Eigen::MatrixXd Sigma = m_Wlist.at(m_depth).selfadjointView<Eigen::Lower>() ;
    return Sigma + std::pow(m_uncorrSD, 2) * Eigen::MatrixXd::Identity(m_Wlist.at(m_depth).rows(), m_Wlist.at(m_depth).cols()) ;
  }
  // The following Kmatrix-related functions work because of the correspondence between knots
  // and observations at the finest resolution.
  mat GetKmatrix() {return m_SigmaInverse ;}
  mat * GetKmatrixAddress() {return &m_SigmaInverse ;}
  mat * GetKmatrixInverseAddress() { return &m_Wlist.at(m_depth) ;}
  mat GetKmatrixInverse() {return GetSigma() - std::pow(m_uncorrSD, 2) * Eigen::MatrixXd::Identity(GetSigma().rows(), GetSigma().cols());}
  vec GetOmega(const uint & order) { throw Rcpp::exception("Trying to get omega vector in tip node! \n") ; return vec(1) ;}
  void SetUncorrSD(const double & sd) {
    m_uncorrSD = sd ;
  }
  mat GetUpred(const uint & l) { return m_UmatList.at(l) ;}
  std::vector<mat> GetUmatList() { return m_UmatList ;}

  void SetPredictLocations(const inputdata &) ;
  Eigen::ArrayXi GetPredIndices() { return m_predsInNode ;}
  void computeUpred(const maternVec &, const maternVec &, const double &, const spatialcoor &, const bool, const double &, const double &) ;

  void genRandomKnots(spatialcoor & dataCoor, const uint & numKnots, const gsl_rng * RNG) {

    m_knotsCoor = spatialcoor(rows(dataCoor.spatialCoords, m_obsInNode),
                              elem(dataCoor.timeCoords, m_obsInNode)) ;
  }

  TipNode(const dimensions & dims, const uint & depth, TreeNode * parent,
          const inputdata & dataset) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_obsInNode = deriveObsInNode(dataset) ;
  }

protected:

  std::vector<mat> GetBlist() {
    return m_Wlist ;
  };

  mat m_SigmaInverse ;
  double m_uncorrSD{ -1 } ;

  // Prediction components (should probably be freed once computations are done)

  // void computeVpred(const vec &, const spatialcoor &) ;
  mat GetLM() { return m_UmatList.at(m_depth) ;}
  // void computeDeltaTildeParas(const inputdata &) ;
  // void recurseBtilde(const uint, const uint) ;
  Eigen::ArrayXi m_predsInNode ;

  // mat m_V ;
  std::vector<mat> m_UmatList ;
  // mat m_covMatrixPredict ;
  // GaussDistParas m_deltaTilde ;
  // std::vector<std::vector<mat>> m_Btilde ;
  // std::vector<mat> m_bPred ;
};
}
