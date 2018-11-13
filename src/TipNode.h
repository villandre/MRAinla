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
  uint GetM() {return m_depth ;}

  void genRandomKnots(inputdata & dataset, uint & numKnots, const gsl_rng * RNG) {
    m_knotsCoor = spatialcoor(dataset.spatialCoords.rows(m_obsInNode),
                              dataset.timeCoords.elem(m_obsInNode)) ;
  }

  TipNode(dimensions & dims, uint & depth, TreeNode * parent, inputdata & dataset, double & covarianceParameter) {
    m_dimensions = dims;
    m_depth = depth ;
    m_parent = parent ;
    m_covPara = covarianceParameter ;
    m_Wlist.resize(m_depth+1) ;
    deriveObsInNode(dataset) ;
  }

protected:
  std::vector<std::vector<arma::mat>> * m_Umat ;
  std::vector<std::vector<arma::mat>> * m_Vmat ;
  std::vector<arma::mat> m_bPred ;
  std::vector<arma::mat> m_bKnots ;
  std::vector<std::vector<arma::mat>> m_Btilde ;
};
}
