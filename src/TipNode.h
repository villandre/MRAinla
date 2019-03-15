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
  uint GetM() {return m_depth ;}

  void DeriveAtilde() {
    for (uint k = 0; k < m_depth; k++) {
      for (uint l = 0; l <= k; l++) {
        m_AtildeList.at(k).at(l) = trans(GetB(k)) * m_SigmaInverse * GetB(l) ;
      }
    }
  }

  void DeriveOmega(const arma::vec &) ;

  void DeriveU(const arma::vec & responseValues) {
    arma::mat uInMat = arma::trans(responseValues.elem(m_obsInNode)) * m_SigmaInverse *
      responseValues.elem(m_obsInNode) ;
    m_u = uInMat.at(0,0) ; // uInMat is supposed to be a 1x1 matrix.
  }

  void DeriveD() {
    double val = 0;
    double sign = 0;
    arma::log_det(val, sign, GetSigma()) ;
    if (sign < 0) {
      throw Rcpp::exception("Sigma is a covariance matrix: it cannot have a negative determinant. \n") ;
    }
    m_d = val ;
  }

  void ComputeWmat(const arma::vec & covParas) {
    baseComputeWmat(covParas) ;
    std::printf("Number of observations: %i \n", m_obsInNode.size()) ;
    std::cout << "Inverting Sigma... \n" ;
    GetSigma().save("/home/luc/Documents/Sigma.info", arma::arma_ascii) ;
    m_SigmaInverse = arma::inv_sympd(GetSigma()) ;
  }

  void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const arma::vec &) ;
  std::vector<std::vector<arma::mat>> GetAlist() const {throw Rcpp::exception("Trying to get an A matrix in a tip node! \n") ;}
  arma::mat GetKtilde() const {throw Rcpp::exception("Trying to get Ktilde in a tip node! \n") ;}
  void deriveBtilde(const spatialcoor & ) ;
  void computeBpred(const spatialcoor &, const arma::vec &) ;
  GaussDistParas CombineEtaDelta(const inputdata &, const arma::vec &) ;
  GaussDistParas GetEtaDelta() const { return m_deltaTilde ;}

  arma::mat GetB(const uint & l) {
    if (l > m_depth) {
      throw Rcpp::exception("Trying to get B^l(j_1, ..., j_M) with l > M! \n") ;
    }
    return m_Wlist.at(l) ; // This works for l = M because knot positions in tips correspond to observation positions. It wouldn't be valid otherwise.
  }

  arma::mat GetSigma() {
    return m_Wlist.at(m_depth) ;
  }
  // The following Kmatrix-related functions work because of the correspondence between knots
  // and observations at the finest resolution.
  arma::mat GetKmatrix() {return m_SigmaInverse ;}
  arma::mat * GetKmatrixAddress() {return &m_SigmaInverse ;}
  arma::mat * GetKmatrixInverseAddress() { return &m_Wlist.at(m_depth) ;}
  arma::mat GetKmatrixInverse() {return GetSigma() ;}
  arma::vec GetOmega(const uint & order) { throw Rcpp::exception("Trying to get omega vector in tip node! \n") ; return arma::vec(1) ;}

  void genRandomKnots(inputdata & dataset, uint & numKnots, const gsl_rng * RNG) {
    m_knotsCoor = spatialcoor(dataset.spatialCoords.rows(m_obsInNode),
                              dataset.timeCoords.elem(m_obsInNode)) ;
  }

  TipNode(const dimensions & dims, const uint & depth, TreeNode * parent,
          const inputdata & dataset) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_obsInNode = deriveObsInNode(dataset) ;
  }

protected:

  std::vector<arma::mat> GetBlist() {
    return m_Wlist ;
  };

  arma::mat m_SigmaInverse ;

  // Prediction components (should probably be freed once computations are done)

  void computeVpred(const arma::vec &, const spatialcoor &) ;
  void computeUpred(const arma::vec &, const spatialcoor &) ;
  arma::mat GetLM() { return m_UmatList.at(m_depth) ;}
  void computeDeltaTildeParas(const inputdata &) ;
  void recurseBtilde(const uint, const uint) ;

  arma::mat m_V ;
  std::vector<arma::mat> m_UmatList ;
  arma::mat m_covMatrixPredict ;
  GaussDistParas m_deltaTilde ;
  std::vector<std::vector<arma::mat>> m_Btilde ;
  std::vector<arma::mat> m_bPred ;
};
}
