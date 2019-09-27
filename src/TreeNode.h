// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::depends(RcppEigen)]]

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_log.h>

#include <RcppGSL.h>
#include <RcppEigen.h>
#include "helper.h"

#ifndef TREENODE_H
#define TREENODE_H

namespace MRAinla
{
typedef unsigned int uint ;
typedef unsigned int uint ;
typedef Eigen::VectorXd vec ;
typedef Eigen::MatrixXd mat ;
typedef Eigen::SparseMatrix<double> sp_mat ;
typedef Eigen::VectorXi uvec ;
typedef Eigen::MatrixXi umat ;

struct maternVec{
  double m_rho ;
  double m_smoothness ;
  double m_scale ;

  friend bool operator==(const maternVec & first, const maternVec & second) {
    return (first.m_rho == second.m_rho) && (first.m_scale == second.m_scale) && (first.m_smoothness == second.m_smoothness) ;
  }

  void print(std::string header) {
    std::cout << header << "\n" ;
    printf("Matern parameters: rho = %.4e smoothness = %.4e scale = %.4e \n", m_rho, m_smoothness, m_scale) ;
  }

  maternVec() {} ;
  maternVec(double rho, double smoothness, double scale) : m_rho(rho), m_smoothness(smoothness), m_scale(scale) { };
};

struct spatialcoor {
  Eigen::ArrayXXd spatialCoords = Eigen::ArrayXXd(1, 2) ;
  Eigen::ArrayXd timeCoords = Eigen::ArrayXd(1) ;

  spatialcoor() { } ;
  spatialcoor(Eigen::ArrayXXd f_sp, Eigen::ArrayXd f_time) : spatialCoords(f_sp), timeCoords(f_time) { } ;

  spatialcoor subset(Eigen::ArrayXi & indices)  const {
    Eigen::ArrayXXd subSpatialCoords = rows(spatialCoords, indices) ;
    Eigen::ArrayXd subTimeCoords = elem(timeCoords, indices) ;
    return spatialcoor(subSpatialCoords, subTimeCoords) ;
  };
};

struct inputdata : public spatialcoor {
  Eigen::VectorXd responseValues = Eigen::VectorXd(1) ;
  Eigen::MatrixXd covariateValues = Eigen::MatrixXd(1, 1) ;

  inputdata() : spatialcoor() {};
  inputdata(Eigen::VectorXd f_responses, Eigen::ArrayXXd f_sp, Eigen::ArrayXd f_time, Eigen::MatrixXd f_covariates)
    : spatialcoor(f_sp, f_time), responseValues(f_responses), covariateValues(f_covariates) {  } ;

  inputdata subset(Eigen::ArrayXi & indices)  const {
    Eigen::ArrayXXd subSpatialCoords = rows(spatialCoords, indices) ;
    Eigen::ArrayXd subTimeCoords = elem(timeCoords, indices) ;
    Eigen::VectorXd subResponseValues = elem(responseValues.array(), indices) ;
    Eigen::MatrixXd subCovariates = rows(covariateValues.array(), indices) ;
    return inputdata(subResponseValues, subSpatialCoords, subTimeCoords, subCovariates) ;
  }
};

struct dimensions {
  Eigen::ArrayXd longitude{Eigen::ArrayXd(2)} ;
  Eigen::ArrayXd latitude{Eigen::ArrayXd(2)} ;
  Eigen::ArrayXd time{Eigen::ArrayXd(2)} ;

  dimensions() { } ;

  dimensions(Eigen::ArrayXd f_lon, Eigen::ArrayXd f_lat, Eigen::ArrayXd f_time)
    : longitude(f_lon), latitude(f_lat), time(f_time) {
    uint firstCompare = (f_lon.size() == f_lat.size()) ;
    uint secondCompare = (f_lat.size() == f_time.size()) ;
    if ((firstCompare * secondCompare) == 0) {
      throw Rcpp::exception("Incompatible data specifications. \n") ;
    }
  } ;
};

struct SpatiotempCoor{
  Eigen::VectorXd sp = Eigen::VectorXd(2) ;
  double time = 0 ;
  SpatiotempCoor(Eigen::VectorXd & sp, double & time) : sp(sp), time(time) { } ;
  SpatiotempCoor() { } ;
};

struct GaussDistParas {
  Eigen::VectorXd meanPara ;
  Eigen::MatrixXd covPara ;
  GaussDistParas() { }
  GaussDistParas(Eigen::VectorXd & meanVec, Eigen::MatrixXd & covMat) : meanPara(meanVec), covPara(covMat) { }
};

class TreeNode
{
public:
  virtual void AddChild(TreeNode *)=0 ;
  virtual void RemoveChild(TreeNode *)=0;
  virtual std::vector<TreeNode *> GetChildren()=0;
  virtual void RemoveChildren()=0;
  virtual int GetM()=0;
  virtual void DeriveAtilde()=0 ;
  virtual void DeriveOmega(const vec &)=0 ;
  virtual void DeriveU(const vec &)=0 ;
  virtual void DeriveD()=0 ;
  virtual void ComputeWmat(const maternVec &, const maternVec &, const double &, const bool, const double &, const double &)=0 ;
  // virtual void ComputeParasEtaDeltaTilde(const spatialcoor &, const inputdata &, const vec &)=0 ;
  virtual std::vector<std::vector<mat>> GetAlist() const = 0;
  virtual mat GetKtilde() const = 0;
  // virtual void deriveBtilde(const spatialcoor & )=0 ;
  // virtual void computeBpred(const spatialcoor &, const vec &)=0 ;
  // virtual GaussDistParas CombineEtaDelta(const inputdata &, const vec &)=0 ;
  // virtual GaussDistParas GetEtaDelta() const =0 ;
  virtual mat GetB(const uint & l)=0 ;
  virtual mat GetSigma()=0 ;
  virtual mat GetKmatrix()=0 ;
  virtual mat * GetKmatrixAddress()=0 ;
  virtual mat * GetKmatrixInverseAddress()=0 ;
  virtual mat GetKmatrixInverse()=0 ;
  virtual vec GetOmega(const uint &)=0 ;
  virtual void SetUncorrSD(const double &)=0 ;
  virtual mat GetUpred(const uint & l)=0 ;
  virtual std::vector<mat> GetUmatList()=0 ;
  virtual void SetPredictLocations(const inputdata &)=0 ;
  virtual Eigen::ArrayXi GetPredIndices()=0 ;
  virtual void computeUpred(const maternVec &, const maternVec &, const double &, const spatialcoor &, const bool, const double &, const double &)=0 ;

  virtual void genRandomKnots(spatialcoor &, const uint &, const gsl_rng *) = 0;

  Eigen::ArrayXi GetObsInNode() {return m_obsInNode ;}
  dimensions GetDimensions() {return m_dimensions;}
  int GetDepth() {return m_depth ;}

  mat GetAtildeList(uint & i, uint & j) {return m_AtildeList.at(i).at(j) ;}
  std::vector<vec> GetOmegaTilde() {return m_omegaTilde ;}
  mat GetOmegaTilde(uint & k) { return m_omegaTilde.at(k) ;}
  spatialcoor GetKnotsCoor() {return m_knotsCoor;}

  std::vector<mat>& GetWlist() {return m_Wlist ;}
  uint GetNodeId() { return m_nodeId ;}

  double GetU() {return m_u ;}
  double GetD() {return m_d ;}
  void SetNodeId(const uint i) { m_nodeId = i ;}

  virtual ~ TreeNode() { } ;

  void SetPredictLocations(const spatialcoor & predictLocations) ;

  Eigen::ArrayXi deriveObsInNode(const spatialcoor &) ;
  void initiateBknots(const vec &) ;
  void completeBknots(const vec &, const uint) ;
  std::vector<mat> GetBknots() const { return m_bKnots ;}


  void clearAtildeList() {m_AtildeList.clear() ;}
  void clearOmegaTilde() {m_omegaTilde.clear() ;}

  uvec GetAncestorIds() {
    std::vector<TreeNode *> ancestorsList = getAncestors() ;
    uvec ancestorIds(ancestorsList.size()) ;
    for (uint i = 0 ; i < ancestorsList.size() ; i++) {
      ancestorIds(i) = ancestorsList.at(i)->GetNodeId() ;
    }
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
  int GetNumKnots() {return m_knotsCoor.timeCoords.size() ;}
  std::vector<TreeNode *> getAncestors() ;

protected:

  TreeNode * m_parent ;
  Eigen::ArrayXi m_obsInNode ;
  int m_depth ;
  dimensions m_dimensions ; // First dimension is longitude, second is latitude, last is time.
  spatialcoor m_knotsCoor ;  // First element is spatial coordinates (longitude, latitude), second is time.
  int m_nodeId ;

  std::vector<std::vector<mat>>& GetAtildeList() {return m_AtildeList ;}
  void baseComputeWmat(const maternVec &, const maternVec &, const double &, const bool, const double &, const double &) ;
  void SetParent(TreeNode * vertexParentPoint) {m_parent = vertexParentPoint ;}

  std::vector<std::vector<mat>> m_AtildeList ;
  std::vector<mat> m_Wlist ;
  std::vector<vec> m_omegaTilde ;
  double m_u ;
  double m_d ;

  double SqExpCovFunction(const Spatiotemprange &, const double &, const double &, const double &, const double &) ;
  double MaternCovFunction(const Spatiotemprange &, const maternVec &, const maternVec &, const double &, const double &, const double &) ;

  mat computeCovMat(const spatialcoor &, const spatialcoor &, const maternVec &, const maternVec &, const double &, const bool, const double &, const double &) ;
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
  // mat ComputeCovMatrix(const vec &) ;

  // For prediction
  void computeBknots() ;

  std::vector<mat> m_bKnots ;
};
}
#endif /* TREENODE_H */
