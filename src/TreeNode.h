// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

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

#include "helper.h"

#ifndef TREENODE_H
#define TREENODE_H

namespace MRAinla
{
typedef unsigned int uint ;

struct spatialcoor {
  arma::mat spatialCoords = arma::mat(1, 2, arma::fill::zeros) ;
  arma::uvec timeCoords = arma::uvec(1, arma::fill::zeros) ;

  spatialcoor() { } ;
  spatialcoor(arma::mat f_sp, arma::uvec f_time) : spatialCoords(f_sp), timeCoords(f_time) { } ;

  spatialcoor subset(arma::uvec & indices)  const {
    arma::mat subSpatialCoords = spatialCoords.rows(indices) ;
    arma::uvec subTimeCoords = timeCoords.rows(indices) ;
    return spatialcoor(subSpatialCoords, subTimeCoords) ;
  };
};

struct inputdata : public spatialcoor {
  arma::vec responseValues = arma::vec(1, arma::fill::zeros) ;
  arma::mat covariateValues = arma::mat(1, 1, arma::fill::zeros) ;

  inputdata() : spatialcoor() {};
  inputdata(arma::vec f_responses, arma::mat f_sp, arma::uvec f_time, arma::mat f_covariates)
    : spatialcoor(f_sp, f_time), responseValues(f_responses), covariateValues(f_covariates) {  } ;

  inputdata subset(arma::uvec & indices)  const {
    arma::mat subSpatialCoords = spatialCoords.rows(indices) ;
    arma::uvec subTimeCoords = timeCoords.elem(indices) ;
    arma::vec subResponseValues = responseValues.elem(indices) ;
    arma::mat subCovariates = covariateValues.rows(indices) ;
    return inputdata(subResponseValues, subSpatialCoords, subTimeCoords, subCovariates) ;
  };
};

struct dimensions {
  arma::vec longitude{arma::vec(2, arma::fill::zeros)} ;
  arma::vec latitude{arma::vec(2, arma::fill::zeros)} ;
  arma::uvec time{arma::uvec(2, arma::fill::zeros)} ;

  dimensions() { } ;

  dimensions(arma::vec f_lon, arma::vec f_lat, arma::uvec f_time)
    : longitude(f_lon), latitude(f_lat), time(f_time) {
    uint firstCompare = (arma::size(f_lon) == arma::size(f_lat)) ;
    uint secondCompare = (arma::size(f_lat) == arma::size(f_time)) ;
    if ((firstCompare * secondCompare) == 0) {
      throw Rcpp::exception("Incompatible data specifications. \n") ;
    }
  } ;
};

struct SpatiotempCoor{
  arma::vec sp = arma::vec(2, arma::fill::zeros) ;
  uint time = 0 ;
  SpatiotempCoor(arma::vec & sp, uint & time) : sp(sp), time(time) { } ;
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
  virtual void DeriveOmega(const inputdata &, const arma::vec &)=0 ;
  virtual void DeriveU(const inputdata &)=0 ;
  virtual void DeriveD()=0 ;
  virtual void ComputeWmat(const arma::vec &)=0 ;
  virtual void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const arma::vec &, const arma::vec &)=0 ;

  virtual void genRandomKnots(inputdata &, uint &, const gsl_rng *) = 0;

  arma::uvec GetObsInNode() {return m_obsInNode ;}
  dimensions GetDimensions() {return m_dimensions;}
  uint GetDepth() {return m_depth ;}

  arma::mat GetAtildeList(uint & i, uint & j) {return m_AtildeList.at(i).at(j) ;}
  arma::mat GetOmegaTilde(uint & k) { return m_omegaTilde.at(k) ;}
  spatialcoor GetKnotsCoor() {return m_knotsCoor;}
  arma::mat GetKmatrix() {return m_K ;}
  arma::mat GetKmatrixInverse() {return m_Kinverse ;}
  std::vector<arma::mat>& GetWlist() {return m_Wlist ;}

  double GetU() {return m_u ;}
  double GetD() {return m_d ;}

  ~ TreeNode() { } ;

  void ComputeBaseKmat(const arma::vec &) ;
  // void SetAtildeList(arma::mat & matrix, uint &i, uint &j) {m_AtildeList.at(i).at(j) = matrix ;}
  void SetPredictLocations(const spatialcoor & predictLocations) ;

  arma::uvec deriveObsInNode(const spatialcoor &) ;

protected:

  TreeNode * m_parent ;
  arma::uvec m_obsInNode ;
  uint m_depth ;
  dimensions m_dimensions ; // First dimension is longitude, second is latitude, last is time.
  spatialcoor m_knotsCoor ;  // First element is spatial coordinates (longitude, latitude), second is time.

  std::vector<std::vector<arma::mat>>& GetAtildeList() {return m_AtildeList ;}
  void baseComputeWmat(const arma::vec &) ;
  TreeNode * GetParent() {return m_parent ;}
  void SetParent(TreeNode * vertexParentPoint) {m_parent = vertexParentPoint ;}

  std::vector<std::vector<arma::mat>> m_AtildeList ;
  std::vector<arma::mat> m_Wlist ;
  std::vector<arma::vec> m_omegaTilde ;
  double m_u ;
  double m_d ;

  arma::mat m_K ;
  arma::mat m_Kinverse ;

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

  // For prediction
  arma::uvec m_predictLocIndices ;
};
}
#endif /* TREENODE_H */
