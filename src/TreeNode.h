// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <RcppEigen.h>

#include "helper.h"

#ifndef TREENODE_H
#define TREENODE_H

namespace MRAinla
{
typedef unsigned int uint ;

struct spatialcoor {
  arma::fmat spatialCoords = arma::fmat(1, 2, arma::fill::zeros) ;
  arma::fvec timeCoords = arma::fvec(1, arma::fill::zeros) ;

  spatialcoor() { } ;
  spatialcoor(arma::fmat f_sp, arma::fvec f_time) : spatialCoords(f_sp), timeCoords(f_time) { } ;

  spatialcoor subset(arma::uvec & indices)  const {
    arma::fmat subSpatialCoords = spatialCoords.rows(indices) ;
    arma::fvec subTimeCoords = timeCoords.rows(indices) ;
    return spatialcoor(subSpatialCoords, subTimeCoords) ;
  };
};

struct inputdata : public spatialcoor {
  arma::vec responseValues = arma::vec(1, arma::fill::zeros) ;
  arma::fmat covariateValues = arma::fmat(1, 1, arma::fill::zeros) ;

  inputdata() : spatialcoor() {};
  inputdata(arma::vec f_responses, arma::fmat f_sp, arma::fvec f_time, arma::fmat f_covariates)
    : spatialcoor(f_sp, f_time), responseValues(f_responses), covariateValues(f_covariates) {  } ;

  inputdata subset(arma::uvec & indices)  const {
    arma::fmat subSpatialCoords = spatialCoords.rows(indices) ;
    arma::fvec subTimeCoords = timeCoords.elem(indices) ;
    arma::vec subResponseValues = responseValues.elem(indices) ;
    arma::fmat subCovariates = covariateValues.rows(indices) ;
    return inputdata(subResponseValues, subSpatialCoords, subTimeCoords, subCovariates) ;
  };
};

struct dimensions {
  arma::fvec longitude{arma::fvec(2, arma::fill::zeros)} ;
  arma::fvec latitude{arma::fvec(2, arma::fill::zeros)} ;
  arma::fvec time{arma::fvec(2, arma::fill::zeros)} ;

  dimensions() { } ;

  dimensions(arma::fvec f_lon, arma::fvec f_lat, arma::fvec f_time)
    : longitude(f_lon), latitude(f_lat), time(f_time) {
    uint firstCompare = (arma::size(f_lon) == arma::size(f_lat)) ;
    uint secondCompare = (arma::size(f_lat) == arma::size(f_time)) ;
    if ((firstCompare * secondCompare) == 0) {
      throw Rcpp::exception("Incompatible data specifications. \n") ;
    }
  } ;
};

struct SpatiotempCoor{
  arma::fvec sp = arma::fvec(2, arma::fill::zeros) ;
  float time = 0 ;
  SpatiotempCoor(arma::fvec & sp, float & time) : sp(sp), time(time) { } ;
  SpatiotempCoor() { } ;
};

struct GaussDistParas {
  arma::vec meanPara ;
  arma::mat covPara ;
  GaussDistParas() { }
  GaussDistParas(arma::vec & meanVec, arma::mat & covMat) : meanPara(meanVec), covPara(covMat) { }
};

class TreeNode
{
public:
  virtual void AddChild(TreeNode *)=0 ;
  virtual void RemoveChild(TreeNode *)=0;
  virtual std::vector<TreeNode *> GetChildren()=0;
  virtual void RemoveChildren()=0;
  virtual uint GetM()=0;
  virtual void DeriveAtilde()=0 ;
  virtual void DeriveOmega(const arma::vec &)=0 ;
  virtual void DeriveU(const arma::vec &)=0 ;
  virtual void DeriveD()=0 ;
  virtual void ComputeWmat(const arma::vec &)=0 ;
  virtual void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const arma::vec &)=0 ;
  virtual std::vector<std::vector<arma::mat>> GetAlist() const = 0;
  virtual arma::mat GetKtilde() const = 0;
  virtual void deriveBtilde(const spatialcoor & )=0 ;
  virtual void computeBpred(const spatialcoor &, const arma::vec &)=0 ;
  virtual GaussDistParas CombineEtaDelta(const inputdata &, const arma::vec &)=0 ;
  virtual GaussDistParas GetEtaDelta() const =0 ;
  virtual arma::mat GetB(const uint & l)=0 ;
  virtual arma::mat GetSigma()=0 ;
  virtual arma::mat GetKmatrix()=0 ;
  virtual arma::mat * GetKmatrixAddress()=0 ;
  virtual arma::mat * GetKmatrixInverseAddress()=0 ;
  virtual arma::mat GetKmatrixInverse()=0 ;

  virtual void genRandomKnots(inputdata &, uint &, const gsl_rng *) = 0;

  arma::uvec GetObsInNode() {return m_obsInNode ;}
  dimensions GetDimensions() {return m_dimensions;}
  uint GetDepth() {return m_depth ;}

  arma::mat GetAtildeList(uint & i, uint & j) {return m_AtildeList.at(i).at(j) ;}
  std::vector<arma::vec> GetOmegaTilde() {return m_omegaTilde ;}
  arma::mat GetOmegaTilde(uint & k) { return m_omegaTilde.at(k) ;}
  spatialcoor GetKnotsCoor() {return m_knotsCoor;}

  std::vector<arma::mat>& GetWlist() {return m_Wlist ;}
  uint GetNodeId() { return m_nodeId ;}

  double GetU() {return m_u ;}
  double GetD() {return m_d ;}
  void SetNodeId(const uint i) { m_nodeId = i ;}

  ~ TreeNode() { } ;

  // void SetAtildeList(arma::mat & matrix, uint &i, uint &j) {m_AtildeList.at(i).at(j) = matrix ;}
  void SetPredictLocations(const spatialcoor & predictLocations) ;

  arma::uvec deriveObsInNode(const spatialcoor &) ;
  void initiateBknots(const arma::vec &) ;
  void completeBknots(const arma::vec &, const uint) ;
  std::vector<arma::mat> GetBknots() const { return m_bKnots ;}

  arma::uvec GetPredictLocIndices() const {return m_predictLocIndices ;}
  void clearAtildeList() {m_AtildeList.clear() ;}
  void clearOmegaTilde() {m_omegaTilde.clear() ;}

  std::vector<uint> GetAncestorIds() {
    std::vector<TreeNode *> ancestorsList = getAncestors() ;
    std::vector<uint> ancestorIds(ancestorsList.size()) ;
    std::transform(ancestorsList.begin(), ancestorsList.end(), ancestorIds.begin(), [] (TreeNode * treeNode) {return treeNode->GetNodeId() ;}) ;
    return ancestorIds ;
  }
  TreeNode * GetParent() {return m_parent ;}

  std::vector<TreeNode *> Siblings() {
    std::vector<TreeNode *> siblingVec ;
    if (m_depth > 0) {
      siblingVec = m_parent->GetChildren() ;
    } else {
      siblingVec.push_back(this) ;
    }
    return siblingVec ;
  }
  uint GetNumKnots() {return m_knotsCoor.timeCoords.size() ;}

protected:

  TreeNode * m_parent ;
  arma::uvec m_obsInNode ;
  uint m_depth ;
  dimensions m_dimensions ; // First dimension is longitude, second is latitude, last is time.
  spatialcoor m_knotsCoor ;  // First element is spatial coordinates (longitude, latitude), second is time.
  uint m_nodeId ;

  std::vector<std::vector<arma::mat>>& GetAtildeList() {return m_AtildeList ;}
  void baseComputeWmat(const arma::vec &) ;
  void SetParent(TreeNode * vertexParentPoint) {m_parent = vertexParentPoint ;}

  std::vector<std::vector<arma::mat>> m_AtildeList ;
  std::vector<arma::mat> m_Wlist ;
  std::vector<arma::vec> m_omegaTilde ;
  double m_u ;
  double m_d ;

  double covFunction(const Spatiotemprange &, const arma::vec &) ;
  std::vector<TreeNode *> getAncestors() ;
  arma::mat computeCovMat(const spatialcoor &, const spatialcoor &, const arma::vec &) ;
  void baseInitialise(const dimensions & dims, const uint & depth, TreeNode * parent, const inputdata & dataset) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    m_Wlist.resize(m_depth+1) ;
    m_AtildeList.resize(m_depth+1) ;
    for (uint i = 0; i < m_AtildeList.size(); i++) {
      m_AtildeList.at(i).resize(i+1) ;
    }
    m_omegaTilde.resize(m_depth + 1) ;
  }
  arma::mat ComputeCovMatrix(const arma::vec &) ;

  // For prediction
  void computeBknots() ;

  arma::uvec m_predictLocIndices ;
  std::vector<arma::mat> m_bKnots ;
};
}
#endif /* TREENODE_H */
