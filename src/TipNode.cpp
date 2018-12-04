#include "TipNode.h"
#include <gsl/gsl_randist.h>

using namespace arma ;
using namespace MRAinla ;

void TipNode::ComputeParasEtaDeltaTilde(const spatialcoor & predictionLocations, const arma::vec & covParas) {
  computeU(predictionLocations, covParas) ;
  computeV(predictionLocations, covParas) ;

}

