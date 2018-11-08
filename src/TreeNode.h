// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]

#include <assert.h>
#include <RcppArmadillo.h>
#include <RcppGSL.h>
#include <gsl/gsl_rng.h>
#include <algorithm>
#include <vector>

using namespace arma ;

#ifndef TREENODE_H
#define TREENODE_H

typedef unsigned int uint ;
typedef std::tuple<vec,vec,uvec> dimtype ;
typedef std::tuple<vec, mat, uvec> datasettype ;
// Is this the ideal way of entering the data? How about std::vector<std::tuple<double, vec, uint>>?
// It would probably make more sense, but the need to compute summaries, e.g. medians, on positions,
// makes this a lot more convenient, even if we don't have a check that the elements are compatible.
typedef std::pair<mat, uvec> knotstype ;

class TreeNode
{
public:
  virtual void AddChild(TreeNode *)=0 ;
  virtual void RemoveChild(TreeNode *)=0;
  virtual std::vector<TreeNode *> GetChildren()=0;
  virtual void RemoveChildren()=0;

  bool IsSolved();
  bool CanSolve();
  void SetSolved(bool);

  void InvalidateSolution();

  std::vector<uvec> * GetInput();

  TreeNode * GetParent() {return _parent ;}
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;}
  void SetId(uint vertexId) {_id = vertexId ;}

  void NegateFlag() {_updateFlag = false ;}

  uint GetId() {return _id ;}
  uvec GetObsInNode() {return _obsInNode ;}
  dimtype GetDimensions() {return _dimensions;}
  uint GetDepth() {return _depth ;}

  ~ TreeNode() { } ;

  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the tree to R.
  TreeNode * _parent ;
  uvec _obsInNode ;
  uint _depth ;
  dimtype _dimensions ; // First dimension is longitude, second is latitude, last is time.
  knotstype _knotsCoor ;  // First element is spatial coordinates (longitude, latitude), second is time.

  void deriveObsInNode(datasettype &) ;

  std::vector<std::vector<mat>> _A ;
  std::vector<std::vector<mat>> _W ;
  std::vector<vec> _omegaTilde ;
  double _u ;
  double _d ;

  mat _K ;
  mat _Kinverse ;

  bool _updateFlag ;
  virtual void genRandomKnots(uint &) = 0;
};

#endif /* TREENODE_H */
