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

class TreeNode
{
public:
  virtual void AddChild(TreeNode *)=0 ;
  virtual void RemoveChild(TreeNode *)=0 ;
  virtual std::vector<TreeNode *> GetChildren()=0;
  virtual void RemoveChildren()=0;

  bool IsSolved();
  bool CanSolve();
  void SetSolved(bool);

  void InvalidateSolution();

  void SetInput(std::vector<uvec> *);
  std::vector<uvec> * GetInput();

  virtual void RestoreIterVecAndExp();

  TreeNode * GetParent() {return _parent ;}
  void SetParent(TreeNode * vertexParentPoint) {_parent = vertexParentPoint ;}
  void SetId(uint vertexId) {_id = vertexId ;}
  uint GetId() {return _id ;}
  void NegateFlag() {_updateFlag = false ;}
  uvec GetObsInNode() {return _obsInNode ;}
  vec GetLonRange() {return _lonRange ;}
  vec GetLatRange() {return _latRange ;}
  uvec GetTimeRange() {return _timeRange ;}

  ~ TreeNode() { } ;

  protected:

  uint _id ; // From 1 to number of nodes. Used for exporting the tree to R.
  TreeNode * _parent ;
  uvec _obsInNode ;
  uint _depth ;
  vec _lonRange ;
  vec _latRange ;
  uvec _timeRange ;
  std::tuple<vec,vec,uvec> dimensions ; // First dimension is longitude, second is latitude, last is time.

  std::vector<std::vector<mat>> _A ;
  std::vector<std::vector<mat>> _W ;
  std::vector<vec> _omegaTilde ;
  double _u ;
  double _d ;

  mat _K ;
  mat _Kinverse ;

  bool _updateFlag ;
};

#endif /* TREENODE_H */
