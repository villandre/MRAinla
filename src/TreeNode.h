// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

#include <algorithm>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include <assert.h>
#include <gsl/gsl_rng.h>

#include <RcppArmadillo.h>
#include <RcppGSL.h>

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
      fputs("Incompatible data specifications \n", stderr) ;
      // exit(EXIT_FAILURE);
    }
  } ;
};

class TreeNode
{
public:
  virtual void AddChild(TreeNode *)=0 ;
  virtual void RemoveChild(TreeNode *)=0;
  virtual std::vector<TreeNode *> GetChildren()=0;
  virtual void RemoveChildren()=0;

  virtual void genRandomKnots(inputdata &, uint &, const gsl_rng *) = 0;

  TreeNode * GetParent() {return m_parent ;}
  void SetParent(TreeNode * vertexParentPoint) {m_parent = vertexParentPoint ;}

  arma::uvec GetObsInNode() {return m_obsInNode ;}
  dimensions GetDimensions() {return m_dimensions;}
  uint GetDepth() {return m_depth ;}
  spatialcoor GetKnotsCoor() {return m_knotsCoor;}

  ~ TreeNode() { } ;

protected:

  TreeNode * m_parent ;
  arma::uvec m_obsInNode ;
  uint m_depth ;
  dimensions m_dimensions ; // First dimension is longitude, second is latitude, last is time.
  spatialcoor m_knotsCoor ;  // First element is spatial coordinates (longitude, latitude), second is time.

  void deriveObsInNode(inputdata &) ;

  std::vector<std::vector<arma::mat>> m_A ;
  std::vector<std::vector<arma::mat>> m_W ;
  std::vector<arma::vec> m_omegaTilde ;
  double m_u ;
  double m_d ;

  arma::mat m_K ;
  arma::mat m_Kinverse ;
};
}
#endif /* TREENODE_H */
