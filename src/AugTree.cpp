#include "AugTree.h"
#include "TipNode.h"
#include "InternalNode.h"
#include "TreeNode.h"

using namespace Rcpp ;
using namespace arma ;

AugTree::AugTree(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, vec & observations, mat & obsSp, uvec & obsTime, uint & minObsForTimeSplit)
{
  _M = M ;
  _dataset = std::make_tuple(observations, obsSp, obsTime) ;
  _mapDimensions = std::make_tuple(lonRange, latRange, timeRange) ;

  BuildTree(minObsForTimeSplit) ;

}

void AugTree::BuildTree(uint & minObsForTimeSplit)
{
  _vertexVector.reserve(1) ;

  // We create the first internal node
  uint depth = 0 ;
  InternalNode * topNode = new InternalNode(_mapDimensions, depth) ;

  createLevels(topNode, minObsForTimeSplit) ;
}

void AugTree::createLevels(TreeNode * parent, uint & numObsForTimeSplit) {

  uvec obsForMedian(std::get<0>(_dataset).size()) ;
  obsForMedian = parent->GetObsInNode() ;
  vec spMedians = median(std::get<1>(_dataset).rows(obsForMedian)) ;
  uint timeMedian = median(std::get<2>(_dataset).elem(obsForMedian)) ;
  std::vector<dimtype> dimensionsForChildren ;
  vec lonRange(2), latRange(2) ;
  uvec timeRange(2) ;

  if (parent->GetObsInNode().size() < numObsForTimeSplit) {
    for (uint i = 0; i < 2; i++) {
      for (uint j = 0 ; j < 2; j++) {
        lonRange.at(0) = std::get<0>(parent->GetDimensions()).at(i) ;
        lonRange.at(1) = spMedians.at(0);
        lonRange = sort(lonRange);

        latRange.at(0) = std::get<1>(parent->GetDimensions()).at(j) ;
        latRange.at(1) = spMedians.at(1);
        latRange = sort(latRange);

        dimensionsForChildren.push_back(std::make_tuple(lonRange, latRange, std::get<2>(parent->GetDimensions()))) ;
      }
    }
  } else {
    for (uint i = 0; i < 2; i++) {
      for (uint j = 0; j < 2; j++) {
        for(uint k = 0; k < 2; k++) {
          lonRange.at(0) = std::get<0>(parent->GetDimensions()).at(i) ;
          lonRange.at(1) = spMedians.at(0);
          lonRange = sort(lonRange);

          latRange.at(0) = std::get<1>(parent->GetDimensions()).at(j) ;
          latRange.at(1) = spMedians.at(1);
          latRange = sort(latRange);

          timeRange.at(0) = std::get<2>(parent->GetDimensions()).at(k) ;
          timeRange.at(1) = timeMedian;
          timeRange = sort(timeRange);

          dimensionsForChildren.push_back(std::make_tuple(lonRange, latRange, timeRange)) ;
        }
      }
    }
  }

  uint incrementedDepth = parent->GetDepth()+1 ;

  for (auto &&i : dimensionsForChildren) {
    if (parent->GetDepth() < _M-1) {
      InternalNode * newNode = new InternalNode(i, incrementedDepth, parent, _dataset) ;
      _vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    } else {
      TipNode * newNode = new TipNode(i, incrementedDepth, parent, _dataset) ;
      _vertexVector.push_back(newNode) ;
      parent->AddChild(newNode) ;
    }
  }

  if (incrementedDepth < _M) {
    incrementedDepth = incrementedDepth + 1 ;
    for (auto && i : parent->GetChildren()) {
      createLevels(i, numObsForTimeSplit) ;
    }
  }
}

void AugTree::InvalidateAll() // Assumes that the tree starts fully solved.
{
  for (auto & i : _vertexVector)
  {
    i->SetSolved(false) ;
  }
}

// void AugTree::BuildEdgeMatrix(uint & M, vec & lonRange, vec & latRange, uvec & timeRange, mat & obsSp, uvec & obsTime, uint & cutForTimeSplit, uint & cutForLonSplit)
// {
//   //TO_DO
// }

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
