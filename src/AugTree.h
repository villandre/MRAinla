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
  dimtype _mapDimensions;

  datasettype _dataset ; // First element is response, second is spatial coordinates, last is time.

  gsl_rng * _randomNumGenerator ;

  // void AddEdgeRecursion(umat &, uint &, TreeNode *) ;
  void BuildTree(uint &) ;
  void createLevels(TreeNode *, uint &) ;
  void generateKnots() ;

public:
  AugTree(uint &, vec &, vec &, uvec &, vec &, mat &, uvec &, uint &) ;

  void InvalidateAll() ;
  // void BuildEdgeMatrix() ;

  void NegateAllUpdateFlags() ;

  std::vector<TreeNode *> GetVertexVector() {return _vertexVector ;} ;

  void ComputeLoglik(const std::vector<mat> &, const std::vector<mat> &, const vec &) ;

  double GetLoglik() {return _logLik ;}
  gsl_rng * GetRandomNumGenerator() {return _randomNumGenerator ;}
  datasettype GetDataset() {return _dataset;}
  knotstype GetKnotsCoor() {return _knotsCoor;}
  uint GetNumTips() {return _numTips ;}

  void InvalidateAllSolutions() ;

  void SetLogLik(double logLik) {_logLik = logLik ;}

  void SetRNG(gsl_rng * myRNG) { _randomNumGenerator = myRNG ;}

  void ComputeLoglik(List &, List &, NumericVector &) ;
  void PrintSolutions(const uint &) ;

  ~AugTree() {deallocate_container(_vertexVector) ;};
};
