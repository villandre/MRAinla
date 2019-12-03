#include "TreeNode.h"

#ifndef MYPROJECT_AUGTREE_H
#define MYPROJECT_AUGTREE_H

struct GammaHyperParas{
  // Default IG mean is 1 (mean = beta/(alpha - 1)), var is beta^2/[(alpha-1)^2(alpha-2)]
  // If alpha <= 2, variance is infinite,
  double m_alpha = 3;
  double m_beta = 2;
  GammaHyperParas() {}
  GammaHyperParas(double alpha, double beta) : m_alpha(alpha), m_beta(beta) { }
  GammaHyperParas(const Eigen::VectorXd & alphaBeta) {
    m_alpha = alphaBeta(0) ;
    m_beta = alphaBeta(1) ;
  }
};

struct maternGammaPriorParasWithoutScale{
  GammaHyperParas m_rho ;
  GammaHyperParas m_smoothness ;

  maternGammaPriorParasWithoutScale() { }
  maternGammaPriorParasWithoutScale(const GammaHyperParas & rho, const GammaHyperParas & smoothness) : m_rho(rho), m_smoothness(smoothness) { }
};

struct maternGammaPriorParas : public maternGammaPriorParasWithoutScale{
  GammaHyperParas m_scale ;

  maternGammaPriorParas(const GammaHyperParas & rho, const GammaHyperParas & smoothness, const GammaHyperParas & scale) : maternGammaPriorParasWithoutScale(rho, smoothness), m_scale(scale) { } ;
};


namespace MRAinla {

class AugTree
{
public:
  AugTree(uint &, uint &, const spaceDimensions &, inputdata &, inputdata &, const unsigned int, double, const std::string &, std::mt19937_64 &) ;

  std::vector<TreeNode *> GetVertexVector() {return m_vertexVector ;} ;

  int GetNumTips() {return m_numTips ;}
  int GetNumKnots() {return m_numKnots ;}
  uint GetNumObs() {return m_numObs ;}

  int GetM() { return m_M ;}

  std::vector<TreeNode *> GetLevelNodes(const uint & level) ;
  std::vector<TreeNode *> GetTipNodes() { return GetLevelNodes(m_M) ;}

  std::vector<mat *> getKmatricesInversePointers() {
    std::vector<mat *> KmatrixInverseList(m_vertexVector.size()) ;
    for (auto & i : m_vertexVector) KmatrixInverseList.at(i->GetNodeId()) = i->GetKmatrixInverseAddress() ;
    return KmatrixInverseList ;
  }
  void computeWmats(const maternVec &, const double &, const double &, const std::string &) ;
  ~ AugTree() {
    deallocate_container(m_vertexVector) ;}

private:

  std::vector<TreeNode *> m_vertexVector ;
  uint m_numObs ;

  int GetNodePos(int nodeId) {
    int nodePos = 0 ;
    while (m_vertexVector.at(nodePos)->GetNodeId() != nodeId) {
      nodePos += 1 ;
    }
  return nodePos ;
  }

  int m_M{ 0 } ;
  int m_Mlon{ 0 } ;
  int m_Mlat{ 0 } ;
  int m_numTips{ 0 } ;
  int m_numKnots{ 0 } ;
  spaceDimensions m_mapDimensions;
  std::string m_distMethod ;

  std::vector<TreeNode *> Descendants(std::vector<TreeNode *>) ;
  void diveAndUpdate(TreeNode *, std::vector<TreeNode *> *) ;

  // Tree construction functions //
  void BuildTree(const unsigned int, double, const inputdata &, const inputdata &, std::mt19937_64 &) ;
  void createLevels(TreeNode *, std::string, Eigen::ArrayXi, const inputdata &) ;
  void generateKnots(TreeNode *, const unsigned int, double, std::mt19937_64 &, const inputdata &, const inputdata &) ;
  void numberNodes() ;
  std::vector<TreeNode *> GetLevel(const uint) ;

  Eigen::Array<bool, Eigen::Dynamic, 1> m_assignedPredToKnot ;
};
}
#endif
