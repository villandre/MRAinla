#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <random>

#include <assert.h>

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
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> sp_mat ;
typedef Eigen::VectorXi uvec ;
typedef Eigen::MatrixXi umat ;
typedef Eigen::Triplet<double> Triplet;

const double epsilon = 1e-15 ;

struct pointerOffset{
  Eigen::MatrixXd * matrixLocation ;
  uint offset ;

  pointerOffset() {
    offset = 0 ;
    matrixLocation = NULL ;
  }
  pointerOffset(Eigen::MatrixXd * matrixLocation, const uint offset):
    matrixLocation(matrixLocation), offset(offset) {} ;
};

struct maternVec{
  double m_rho ;
  double m_smoothness ;
  double m_scale ;

  friend bool operator==(const maternVec & first, const maternVec & second) {
    return (fabs(first.m_rho - second.m_rho) < epsilon) && (fabs(first.m_scale - second.m_scale) < epsilon) && (fabs(first.m_smoothness - second.m_smoothness) < epsilon) ;
  }

  void print(std::string header) const {
    Rcpp::Rcout << header << "\n" ;
    Rprintf("Matern parameters: rho = %.4e smoothness = %.4e scale = %.4e \n", m_rho, m_smoothness, m_scale) ;
  }

  maternVec() {
    m_rho = 0;
    m_smoothness = 0 ;
    m_scale = 0 ;
  } ;
  maternVec(double rho, double smoothness, double scale) : m_rho(rho), m_smoothness(smoothness), m_scale(scale) { };
};

struct spatialcoor {
  Eigen::ArrayXXd spatialCoords  ;

  spatialcoor() {
    spatialCoords = Eigen::ArrayXXd(0, 2) ;
  } ;
  spatialcoor(Eigen::ArrayXXd f_sp) : spatialCoords(f_sp) { } ;

  spatialcoor subset(Eigen::ArrayXi & indices)  const {
    Eigen::ArrayXXd subSpatialCoords = rows(spatialCoords, indices) ;
    return spatialcoor(subSpatialCoords) ;
  };
  void print() {
    Rcpp::Rcout << "Space coordinates:" << std::endl ;
    Rcpp::Rcout << spatialCoords << std::endl ;
  }
};

struct spatiotempcoor : public spatialcoor {
  Eigen::ArrayXd timeCoords ;
  spatiotempcoor() : spatialcoor() { } ;
  spatiotempcoor(Eigen::ArrayXXd spatialCoords, Eigen::ArrayXd timeCoords) : spatialcoor(spatialCoords),
  timeCoords(timeCoords) {  }
  void print() {
    Eigen::ArrayXXd merged(spatialCoords.rows(), 3) ;
    merged << spatialCoords, timeCoords ;
    Rcpp::Rcout << "Spatiotemporal coords: " << std::endl ;
    Rcpp::Rcout << merged << std::endl ;
  }
};

struct inputdata : public spatiotempcoor {
  Eigen::VectorXd responseValues ;
  Eigen::MatrixXd covariateValues ;

  inputdata() : spatiotempcoor() {
    responseValues = Eigen::VectorXd(0) ;
    covariateValues = Eigen::MatrixXd(0, 0) ;
  };

  inputdata(Eigen::VectorXd f_responses, Eigen::ArrayXXd f_sp, Eigen::ArrayXd f_time, Eigen::MatrixXd f_covariates)
    : spatiotempcoor(f_sp, f_time), responseValues(f_responses), covariateValues(f_covariates) {  } ;

  inputdata subset(Eigen::ArrayXi & indices)  const {
    Eigen::ArrayXXd subSpatialCoords = rows(spatialCoords, indices) ;
    Eigen::ArrayXd subTimeCoords = elem(timeCoords, indices) ;
    Eigen::VectorXd subResponseValues = elem(responseValues.array(), indices) ;
    Eigen::MatrixXd subCovariates = rows(covariateValues.array(), indices) ;
    return inputdata(subResponseValues, subSpatialCoords, subTimeCoords, subCovariates) ;
  }
};

struct spaceDimensions {
  Eigen::Array2d longitude ;
  Eigen::Array2d latitude ;

  void print() {
    Rcpp::Rcout << "Longitude range: \n" << longitude << "\n\n" ;
    Rcpp::Rcout << "Latitude range: \n" << latitude << "\n\n" ;
  }

  spaceDimensions() {
    longitude = Eigen::Array2d::Zero(2) ;
    latitude = Eigen::Array2d::Zero(2) ;
  } ;

  spaceDimensions(Eigen::ArrayXd f_lon, Eigen::ArrayXd f_lat)
    : longitude(f_lon), latitude(f_lat) {
    uint firstCompare = (f_lon.size() == f_lat.size()) ;
    if (firstCompare == 0) {
      throw Rcpp::exception("Incompatible data specifications. \n") ;
    }
  } ;
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
  virtual void ComputeWmat(const maternVec &, const double &, const double &, const std::string &)=0 ;
  virtual mat & GetB(const uint & l)=0 ;
  virtual mat GetSigma()=0 ;
  virtual mat & GetKmatrix()=0 ;
  virtual mat * GetKmatrixAddress()=0 ;
  virtual mat * GetKmatrixInverseAddress()=0 ;
  virtual mat GetKmatrixInverse()=0 ;
  virtual void SetUncorrSD(const double &)=0 ;
  virtual mat & GetUpred(const uint & l)=0 ;
  virtual std::vector<mat> & GetUmatList()=0 ;
  virtual void SetPredictLocations(const inputdata &)=0 ;
  virtual Eigen::ArrayXi & GetPredIndices()=0 ;
  virtual void computeUpred(const maternVec &, const double &, const spatialcoor &, const double &, const std::string &)=0 ;

  virtual void genKnotsOnSquare(spatialcoor &, int &, std::mt19937_64 &, Eigen::Array<bool, Eigen::Dynamic, 1> &) = 0;
  virtual void genRandomKnots(spatialcoor &, int &, std::mt19937_64 &) = 0;

  void clearWmatrices() {
    for (auto & i : m_Wlist) {
      i.resize(0, 0) ;
    }
  }

  Eigen::ArrayXi & GetObsInNode() {return m_obsInNode ;}
  spaceDimensions GetDimensions() {return m_dimensions;}
  int GetDepth() {return m_depth ;}

  spatialcoor & GetKnotsCoor() {return m_knotsCoor;}

  std::vector<mat> & GetWlist() {return m_Wlist ;}
  uint GetNodeId() { return m_nodeId ;}
  void SetNodeId(const uint i) { m_nodeId = i ;}

  virtual ~ TreeNode() { } ;

  void SetPredictLocations(const spatialcoor & predictLocations) ;

  Eigen::ArrayXi deriveObsInNode(const spatialcoor &) ;

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
  int GetNumKnots() {return m_knotsCoor.spatialCoords.rows();}
  std::vector<TreeNode *> getAncestors() ;
  int GetNumObs() {return m_obsInNode.size() ;}

protected:

  TreeNode * m_parent ;
  Eigen::ArrayXi m_obsInNode ;
  int m_depth{ -1 } ;
  spaceDimensions m_dimensions ; // First dimension is longitude, second is latitude
  spatialcoor m_knotsCoor ;
  int m_nodeId{ -1 } ;

  void baseComputeWmat(const maternVec &, const double &, const double &, const std::string &) ;
  void SetParent(TreeNode * vertexParentPoint) {m_parent = vertexParentPoint ;}

  std::vector<mat> m_Wlist ;

  double MaternCovFunction(const double &, const maternVec &, const double &, const double &) ;

  mat computeCovMat(const spatialcoor &, const spatialcoor &, const maternVec &, const double &, const double &, const std::string &) ;
  void baseInitialise(const spaceDimensions & dims, const uint & depth, TreeNode * parent, const inputdata & dataset) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    m_Wlist.resize(m_depth+1) ;
  }

  // For prediction
  void computeBknots() ;

  std::vector<mat> m_bKnots ;
};
}
#endif /* TREENODE_H */
