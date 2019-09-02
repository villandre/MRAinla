#include "TreeNode.h"

#ifndef INTERMEDIATENODE_H
#define INTERMEDIATENODE_H
namespace MRAinla
{
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
  void DeriveAtilde() ;
  void DeriveOmega(const arma::vec &) ;
  void DeriveU(const arma::vec &) ;
  void DeriveD() ;
  void ComputeWmat(const maternVec &, const maternVec &, const double &, const bool, const double &, const double &) ;
  void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const arma::vec &) ;
  std::vector<std::vector<arma::mat>> GetAlist() const {return m_Alist ;};
  arma::mat GetKtilde() const {return m_Ktilde ;}
  // void deriveBtilde(const spatialcoor & x1) { throw Rcpp::exception("Btilde are only needed for tips. \n") ;}
  // void computeBpred(const spatialcoor & x1, const arma::vec & x2) { throw Rcpp::exception("Bpred should only be computed in tips. \n") ;}
  // GaussDistParas CombineEtaDelta(const inputdata &, const arma::vec &) { throw Rcpp::exception("Combination should only occur in tip nodes. \n") ;}
  // GaussDistParas GetEtaDelta() const { return m_etaTilde ;}
  arma::mat GetB(const uint & l) { throw Rcpp::exception("Trying to get B matrix in internal node.\n") ;}
  arma::mat GetSigma() { throw Rcpp::exception("Trying to get Sigma matrix in internal node.\n") ;}
  arma::mat GetKmatrix() {return m_K ;}
  arma::mat * GetKmatrixAddress() {return &m_K ;}
  arma::mat * GetKmatrixInverseAddress() { return &(m_Wlist.back()) ;}
  arma::mat GetKmatrixInverse() {return m_Wlist.back() ;}
  arma::vec GetOmega(const uint & order) { return m_omega.at(order) ;}
  void SetUncorrSD(const double &) {throw Rcpp::exception("Trying to add uncorrelated error for internal nodes! \n") ;}
  arma::mat GetUpred(const uint & l) { throw Rcpp::exception("Upred matrices only computed for tip nodes! \n") ;}
  std::vector<arma::mat> GetUmatList() { throw Rcpp::exception("UmatList only available in tip nodes! \n") ;}
  void SetPredictLocations(const inputdata & data) { throw Rcpp::exception("Trying to attach predict locations to internal nodes! Predict locations should only be defined in the tips! \n") ;}
  arma::uvec GetPredIndices() { throw Rcpp::exception("Prediction locations not defined in internal nodes! \n");}
  void computeUpred(const maternVec &, const maternVec &, const double &, const spatialcoor &, const bool, const double &, const double &) {
    throw Rcpp::exception("Upred matrices need not be computed in internal nodes! \n") ;
  }

  void genRandomKnots(spatialcoor &, const uint &, const uint &, const gsl_rng *) ;

  InternalNode(const dimensions & dims, const uint & depth, TreeNode * parent,
               const inputdata & dataset) {
    baseInitialise(dims, depth, parent, dataset) ;
    m_obsInNode = deriveObsInNode(dataset) ;
    m_omega.resize(m_depth + 1) ;
    m_Alist.resize(m_depth+1) ;
    for (uint i = 0; i < m_Alist.size(); i++) {
      m_Alist.at(i).resize(i+1) ;
    }
  }

  InternalNode(const dimensions & dims, const inputdata & dataset) {
    baseInitialise(dims, 0, this, dataset) ;
    int numObs = dataset.responseValues.size() ;
    m_obsInNode = arma::regspace<arma::uvec>(0, numObs - 1) ;
    m_omega.resize(m_depth + 1) ;
    m_Alist.resize(m_depth+1) ;
    for (uint i = 0; i < m_Alist.size(); i++) {
      m_Alist.at(i).resize(i+1) ;
    }
  }

protected:

  std::vector<TreeNode *> m_children ;
  std::vector<std::vector<arma::mat>> m_Alist ;
  std::vector<arma::vec> m_omega ;
  arma::mat m_Ktilde ;
  arma::mat m_KtildeInverse ;
  arma::mat m_K ;

  // Prediction elements

  GaussDistParas m_etaTilde ;
};
}
#endif /* INTERMEDIATENODE_H */
