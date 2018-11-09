#include "TreeNode.h"

namespace MRAinla
{
class TipNode:public TreeNode
{
public:

  void AddChild(TreeNode * child) {assert(false) ;}
  void RemoveChild(TreeNode *) {assert(false) ;}
  std::vector<TreeNode *> GetChildren() {
    std::vector<TreeNode *> myVec;
    myVec.push_back(NULL) ;
    return myVec;
  } // An input node returns a null pointer when it is asked to provide the address of a child.
  void RemoveChildren() {}

  bool IsSolved() {return true ;}
  bool CanSolve() {return true ;}
  void SetSolved(bool status) {}

  TipNode(dimtype & dims, uint & depth, TreeNode * parent, datasettype & dataset) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    deriveObsInNode(dataset) ;
  }


protected:
  std::vector<std::vector<arma::mat>> * m_Umat ;
  std::vector<std::vector<arma::mat>> * m_Vmat ;
  std::vector<arma::mat> m_bPred ;
  std::vector<arma::mat> m_bKnots ;
  std::vector<std::vector<arma::mat>> m_Btilde ;
  void genRandomKnots(datasettype & dataset, uint & numKnots, gsl_rng * RNG) {
    m_knotsCoor = std::make_tuple(std::get<1>(dataset).rows(m_obsInNode),
                                  std::get<2>(dataset).elem(m_obsInNode)) ;
  }
};
}
