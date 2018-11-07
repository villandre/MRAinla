#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"
#include "TreeNode.h"

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, vec & observations, mat & obsSp, uvec & obsTime, uint & minObsForLonSplit, uint & minObsForTimeSplit)
{
  _M = M ;
  _responseValues = observations ;
  _obsSp = obsSp ;
  _obsTime = obsTime ;
  _lonRange = lonRange ;
  _latRange = latRange ;
  _timeRange = timeRange ;

  BuildTree(minObsForLonSplit, minObsForTimeSplit) ;

}

void AugTree::BuildTree(const uint & minObsForLonSplit, const uint & minObsForTimeSplit)
{
  _vertexVector.reserve(1) ;

  // We create the first internal node

  InternalNode * topNode = new InternalNode(_lonRange, _latRange, _timeRange, _obsSp, _obsTime, 0) ;
  topNode->SetParent(topNode) ; // The root node has itself as a parent



  createLevels(1, topNode) ;
  for (uint level = 1; level <= _M; level++) {


    splitSizes

    if (level == _M) {
      TipNode * newNode = new TipNode() ;
    } else {
      IntermediateNode *
    }
    InternalNode * newNode = new TipNode() ;
    _vertexVector.push_back(newNode) ;
  }

  // We add the internal nodes
  for (uint i = 0 ; i < edgeMatrix.n_rows - _numTips + 1; i++) {
    InternalNode * newNode = new InternalNode() ;
    _vertexVector.push_back(newNode) ;
  }
  // We set the IDs (to facilitate exporting the phylogeny to R).
  for (uint i = 0 ; i < _vertexVector.size(); i++) {
    _vertexVector[i]->SetId(i) ;
  }
  // The vertices are all disjoint, the next loop defines their relationships
  // The iterator follows columns.
  for (umat::const_iterator iter = edgeMatrix.begin(); iter < edgeMatrix.end()-edgeMatrix.n_rows; iter++)
  {
    _vertexVector[*iter]->AddChild(_vertexVector[*(iter+edgeMatrix.n_rows)]) ;
    _vertexVector[*(iter+edgeMatrix.n_rows)]->SetParent(_vertexVector[*iter]) ;
  }
}

void AugTree::createLevels(uint & depth, TreeNode * parent) {

  uvec obsForMedian(_responseValues.size()) ;
  obsForMedian = parent->GetObsInNode() ;
  vec spMedians = median(_obsSp.rows(obsForMedian)) ;
  uint timeMedian = median(_obsTime.elem(obsForMedian)) ;
  std::tuple<vec,vec,uvec> dimensions(8) ;

  for (uint i = 0 ; i < 8; i++) {

    std::tuple<vec,vec,uvec> components(1) ;
    vec longitudeRange(2), latitudeRange(2) ;
    uvec timeRange(2) ;
    size_t index ;

    std::bitset<3> combination(i) ;

    longitudeRange.at(0) = spMedians.at(0) ;
    index = (combination[0]) ; // Typecasting to integer 0 or 1.
    longitudeRange.at(1) = parent->GetLonRange().at(index) ;
    longitudeRange = sort(longitudeRange) ;

    latitudeRange.at(0) = spMedians.at(1) ;
    index = (combination[1]) ;
    latitudeRange.at(1) = parent->GetLatRange().at(index) ;
    latitudeRange = sort(latitudeRange) ;

    timeRange.at(0) = timeMedian ;
    index = (combination[1]) ;
    latitudeRange.at(1) = parent->GetLatRange().at(index) ;
    latitudeRange = sort(latitudeRange) ;
  }



  if (depth < _M) {
    createLevel(depth+1, )
  }
}

void AugTree::InvalidateAll() // Assumes that the tree starts fully solved.
{
  for (auto & i : _vertexVector)
  {
    i->SetSolved(false) ;
  }
}

void AugTree::BuildEdgeMatrix(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, mat & obsSp, uvec & obsTime, uint & cutForTimeSplit, uint & cutForLonSplit)
{
  for (uint i = 0; i == M; i++) {
    double lonMedian = median(obsSp.col(0)) ;
    double latMedian = median(obsSp.col(1)) ;
    uint timeMedian = median(obsTime) ;

  }


}

void AugTree::AddEdgeRecursion(umat & matToUpdate, uint & lineNum, TreeNode * currentVertex)
{
  if (currentVertex->GetParent() != NULL)
  {
    matToUpdate.at(lineNum, 0) = currentVertex->GetParent()->GetId()+1 ; // Vertex numbering starts at 1 in R, instead of 0.
    matToUpdate.at(lineNum, 1) = currentVertex->GetId()+1 ;
    lineNum = lineNum + 1 ;
  }
  if (currentVertex->GetChildren().at(0) != NULL)
  {
    for (auto & i : currentVertex->GetChildren())
    {
      AddEdgeRecursion(matToUpdate, lineNum, i) ;
    }
  }
}

void AugTree::ComputeLoglik(const std::vector<mat> & withinClusTransProbs, const std::vector<mat> & betweenClusTransProbs, const vec & limProbs)
{
// TO_DO
}

void AugTree::NegateAllUpdateFlags()
{
  for (auto & i : _vertexVector)
  {
    i->NegateFlag() ;
  }
}

void AugTree::PrintSolutions(const uint & elementNum)
{
  // for (auto & i : _vertexVector)
  // {
  //   cout << "This is node " << i->GetId() << ".";
  //   cout << "Is my supporting branch within-cluster? " << i->GetWithinParentBranch() << ".\n";
  //   i->GetDictionaryIterator(elementNum, _numRateCats)->second.first.print("Solution from dictionary (could be rescaled):") ;
  // }
}
