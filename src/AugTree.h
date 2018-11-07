#include "TreeNode.h"
#include "helper.h"

using namespace arma ;
using namespace Rcpp ;

class AugTree
{
protected:
  double _logLik ;
  std::vector<TreeNode *> _vertexVector ;

  uint _M ;
  uint _numTips ;
  umat _edgeMatrix ;
  vec _lonRange ;
  vec _latRange ;
  uvec _timeRange ;

  vec _responseValues ;
  mat _obsSp ;
  uvec _obsTime ;

  std::vector<vec> _knotsSp ;
  DateVector _knotsTime ;

  gsl_rng * _randomNumGenerator ;

  void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  void BuildTree(const uint &, const uint &) ;
  void createLevels(uint &, TreeNode *) ;

public:
  AugTree(uint &, vec &, vec &, uvec &, vec &, mat &, uvec &, uint &, uint &) ;

  void InvalidateAll() ;
  void BuildEdgeMatrix() ;

  void NegateAllUpdateFlags() ;

  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;

  void ComputeLoglik(const std::vector<mat> &, const std::vector<mat> &, const vec &) ;
  double GetLoglik() {return _logLik ;}
  gsl_rng * GetRandomNumGenerator() {return _randomNumGenerator ;}

  uint GetNumTips() {return _numTips ;}

  void InvalidateAllSolutions() ;

  void SetLogLik(double logLik) {_logLik = logLik ;}

  void SetRNG(gsl_rng * myRNG) { _randomNumGenerator = myRNG ;}

  void ComputeLoglik(List &, List &, NumericVector &) ;
  void PrintSolutions(const uint &) ;

  ~AugTree() {deallocate_container(_vertexVector) ;};
};
