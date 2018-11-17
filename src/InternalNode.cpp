#include "InternalNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void InternalNode::RemoveChild(TreeNode * childToRemove)
{
  auto ChildIterPos = std::find(m_children.begin(), m_children.end(), childToRemove) ;
  if (ChildIterPos == m_children.end())
  {
    cerr << "Warning: Trying to remove a child that was not found! \n" ;
  }
  else
  {
    m_children.erase(ChildIterPos) ;
  }
}

void InternalNode::genRandomKnots(inputdata & dataset, uint & numKnots, const gsl_rng * RNG) {

  mat knotsSp(numKnots, 2, fill::zeros) ;

  double minLon = min(m_dimensions.longitude) ;
  double maxLon = max(m_dimensions.longitude) ;

  double minLat = min(m_dimensions.latitude) ;
  double maxLat = max(m_dimensions.latitude) ;

  for (mat::iterator iter = knotsSp.begin() ; iter != (knotsSp.begin() + knotsSp.n_rows) ; iter++) {
    (*iter) = gsl_ran_flat(RNG, minLon, maxLon) ;
    *(iter + knotsSp.n_rows) = gsl_ran_flat(RNG, minLat, maxLat) ;
  }

  double minTime = (double) min(m_dimensions.time) ;
  uint timeRangeSize = range(m_dimensions.time) ;
  uvec time(numKnots) ;

  time.imbue( [&]() { return gsl_rng_uniform_int(RNG, timeRangeSize) + minTime; } ) ;
  m_knotsCoor = spatialcoor(knotsSp, time) ;
}

void InternalNode::DeriveAtilde() {
  cout << "Entering DeriveAtilde... Entering loop... \n" ;
  mat containerMat ;
  for (uint k = 0; k <= m_depth ; k++) {
    for (uint l = 0; l <= k ; l++) {
      containerMat.reshape(m_children.at(0)->GetAtildeList(k, l).n_rows,
                           m_children.at(0)->GetAtildeList(k, l).n_cols) ;
      containerMat.fill(0) ;

      // containerMat = std::accumulate(m_children.begin(), m_children.end(), containerMat,
      //                                [&k, &l](const mat & a, TreeNode * b) {
      //                                  return a + b->GetAtildeList(k, l);
      //                                }) ;
      for (auto & i : m_children) {
        containerMat += i->GetAtildeList(k,l) ;
      }
      m_Alist.at(k).at(l) = containerMat ;
    }
  }
  cout << "Cleared loop! \n Computing KtildeInverse..." ;
  mat KtildeInverse = m_Kinverse + m_Alist.at(m_depth).at(m_depth) ;
  m_KtildeInverse = KtildeInverse ;
  mat Ktilde = symmatl(KtildeInverse) ;
  m_Ktilde = Ktilde ;
  cout << "Entering second loop... \n" ;
  for (uint k = 0; k <= m_depth ; k++) {
    for (uint l = 0; l <= k ; l++) {
      // printf("First term size: %u, %u \n", m_Alist.at(k).at(l).n_rows, m_Alist.at(k).at(l).n_cols) ;
      // printf("Second term first mat size: %u, %u \n", m_Alist.at(k).at(m_depth).n_rows, m_Alist.at(k).at(m_depth).n_cols) ;
      // printf("Second term middle mat size: %u, %u \n", m_Ktilde.n_rows, m_Ktilde.n_cols) ;
      // printf("Second term last mat size: %u, %u \n", m_Alist.at(m_depth).at(l).n_rows, m_Alist.at(m_depth).at(l).n_cols) ;
      mat Atilde = m_Alist.at(k).at(l) - m_Alist.at(k).at(m_depth) * m_Ktilde *
        m_Alist.at(m_depth).at(l) ;
      m_AtildeList.at(k).at(l) = Atilde ;
    }
  }
}

void InternalNode::DeriveOmega(const inputdata & dataset) {

  vec containerVec ;
  for (uint k = 0; k <= m_depth ; k++) {
    containerVec.resize(m_children.at(0)->GetOmegaTilde(k).size()) ;
    containerVec.fill(0) ;
    containerVec = std::accumulate(m_children.begin(), m_children.end(), containerVec,
                                   [&k](vec a, TreeNode * b) {
                                     return a + b->GetOmegaTilde(k);
                                   }) ;
    m_omega.at(k) = containerVec ;
    m_omegaTilde.at(k) = m_omega.at(k) - m_Alist.at(k).at(m_depth) * m_Ktilde * m_omega.at(m_depth) ;
  }
}

void InternalNode::DeriveU(const inputdata & dataset) {
  mat firstTerm = -trans(m_omega.at(m_depth)) * m_Ktilde * m_omega.at(m_depth) ;
  double secondTerm = 0 ;
  secondTerm = std::accumulate(m_children.begin(), m_children.end(), secondTerm,
                               [](double a, TreeNode * b) {
                                 return a + b->GetU();
                               }) ;
  m_u = firstTerm.at(0,0) + secondTerm ;
}

void::InternalNode::DeriveD() {
  double value = 0 ;
  double sign = 0 ;
  log_det(value, sign, m_KtildeInverse) ; // Function is supposed to update value and sign in place.
  m_d = gsl_sf_log(sign) + value ;
  log_det(value, sign, m_Kinverse) ;
  m_d = m_d - (gsl_sf_log(sign) + value) ;
  double thirdTerm = 0 ;
  thirdTerm = std::accumulate(m_children.begin(), m_children.end(), thirdTerm,
                               [](double a, TreeNode * b) {
                                 return a + b->GetD();
                               }) ;
 m_d = m_d + thirdTerm ;
}

