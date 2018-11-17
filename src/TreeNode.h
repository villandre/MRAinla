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
};

struct inputdata : public spatialcoor {
  arma::vec responseValues = arma::vec(1, arma::fill::zeros) ;

  inputdata() : spatialcoor() {};
  inputdata(arma::vec f_responses, arma::mat f_sp, arma::uvec f_time)
    : spatialcoor(f_sp, f_time), responseValues(f_responses) {  } ;
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

class TreeNode
{
public:
  virtual void AddChild(TreeNode *)=0 ;
  virtual void RemoveChild(TreeNode *)=0;
  virtual std::vector<TreeNode *> GetChildren()=0;
  virtual void RemoveChildren()=0;
  virtual uint GetM()=0;
  virtual void DeriveAtilde()=0 ;
  virtual void DeriveOmega(const inputdata &)=0 ;
  virtual void DeriveU(const inputdata &)=0 ;
  virtual void DeriveD()=0 ;
  virtual void DeriveB()=0 ;

  virtual void genRandomKnots(inputdata &, uint &, const gsl_rng *) = 0;

  TreeNode * GetParent() {return m_parent ;}
  void SetParent(TreeNode * vertexParentPoint) {m_parent = vertexParentPoint ;}

  arma::uvec GetObsInNode() {return m_obsInNode ;}
  dimensions GetDimensions() {return m_dimensions;}
  uint GetDepth() {return m_depth ;}
  spatialcoor GetKnotsCoor() {return m_knotsCoor;}
  arma::mat GetKmatrix() {return m_K ;}
  arma::mat GetKmatrixInverse() {return m_Kinverse ;}
  std::vector<arma::mat> GetWlist() {return m_Wlist ;}
  std::vector<arma::mat> GetBlist() {return m_Blist ;}
  arma::mat GetAtildeList(uint & i, uint & j) {return m_AtildeList.at(i).at(j) ;}
  std::vector<std::vector<arma::mat>> GetAtildeList() {return m_AtildeList ;}
  arma::mat GetOmegaTilde(uint & k) { return m_omegaTilde.at(k) ;}
  double GetU() {return m_u ;}
  double GetD() {return m_d ;}

  ~ TreeNode() { } ;
  void ComputeWmat() ;
  void ComputeBaseKmat() ;
  void SetSigma(arma::mat & SigmaMatrix) {
    m_Sigma = SigmaMatrix ;
    m_SigmaInverse = inv_sympd(SigmaMatrix) ;
  }
  void SetAtildeList(arma::mat & matrix, uint &i, uint &j) {m_AtildeList.at(i).at(j) = matrix ;}

protected:

  TreeNode * m_parent ;
  arma::uvec m_obsInNode ;
  uint m_depth ;
  dimensions m_dimensions ; // First dimension is longitude, second is latitude, last is time.
  spatialcoor m_knotsCoor ;  // First element is spatial coordinates (longitude, latitude), second is time.

  void deriveObsInNode(const inputdata &) ;

  std::vector<std::vector<arma::mat>> m_AtildeList ;
  std::vector<arma::mat> m_Wlist ;
  std::vector<arma::vec> m_omegaTilde ;
  std::vector<arma::mat> m_Blist ;
  double m_u ;
  double m_d ;

  arma::mat m_K ;
  arma::mat m_Kinverse ;
  arma::mat m_Sigma ;
  arma::mat m_SigmaInverse ;

  arma::vec m_covPara ;

  double covFunction(const Spatiotemprange &) ;
  std::vector<TreeNode *> getAncestors() ;
  arma::mat computeCovMat(const spatialcoor &, const spatialcoor &) ;
  void baseInitialise(const dimensions & dims, const uint & depth, TreeNode * parent, const inputdata & dataset, const arma::vec & covPars) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    m_covPara = covPars ;
    m_Wlist.resize(m_depth+1) ;
    m_AtildeList.resize(m_depth+1) ;
    for (uint i = 0; i < m_AtildeList.size(); i++) {
      m_AtildeList.at(i).resize(i+1) ;
    }
    m_Blist.resize(m_depth + 1) ;
    m_omegaTilde.resize(m_depth + 1) ;
  }
};
}
#endif /* TREENODE_H */
