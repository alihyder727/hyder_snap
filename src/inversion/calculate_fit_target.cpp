/** @file calculate_fit_target.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Nov 18, 2021 21:07:59 EST
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>
#include <iomanip>

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../radiation/radiation.hpp"
#include "../utils/utils.hpp"
#include "../debugger/debugger.hpp"
#include "../math/linalg.h"

void calculate_fit_target(MeshBlock *pmb, Real *val, int nvalue,
    int k, int j, bool differential)
{
  //ATHENA_LOG("calculate_fit_target");
  std::stringstream &msg = pmb->pdebug->msg;
  std::vector<Direction> out_dir = pmb->prad->GetOutgoingRays();
  int nangle = 0;

  pmb->pdebug->Call("calculate_fit_target");
  msg << "- model " << j - pmb->js << std::endl;
  msg << "- emission angles used: ";
  for (int i = 0; i < out_dir.size(); ++i)
    if (out_dir[i].mu >= cos(45./180.*M_PI)) {
      msg << acos(out_dir[i].mu)/M_PI*180 << " ";
      nangle++;
    }
  msg << std::endl;

  Real **B, *b;
  NewCArray(B, nangle, 3);
  b = new Real [nangle];
  
  // 10. solve coefficients using least square
  for (int i = 0; i < nangle; ++i) {
    B[i][0] = 1.;
    B[i][1] = 1. - out_dir[i].mu;
    B[i][2] = (1. - out_dir[i].mu)*(1. - out_dir[i].mu);
  }

  // 11. log likelihood
  RadiationBand *pband = pmb->prad->pband;
  int i = 0;
  while (pband != NULL) {
    for (int n = 0; n < nangle; ++n)
      b[n] = pband->btoa(n,k,j);
    leastsq(B, b, nangle, 3);
    for (int n = 0; n < 3; ++n)
      val[i*3+n] = b[n];

    if (differential) {
      for (int n = 0; n < nangle; ++n)
        b[n] = pband->btoa(n,pmb->ks,pmb->js);
      leastsq(B, b, nangle, 3);
      for (int n = 0; n < 3; ++n)
        val[i*3+n] -= b[n];
    }

    i++;
    pband = pband->next;
  }

  if (nvalue != 3*i) {
    msg << "### FATAL ERROR in calculate_fit_target" << std::endl
        << "Lengths of forward vector and target vector do not match";
    ATHENA_ERROR(msg);
  }

  msg << "- foward model results: ";
  for (int i = 0; i < nvalue; ++i)
    msg << std::setprecision(5) << val[i] << " ";
  msg << std::endl;

  FreeCArray(B);
  delete[] b;
  pmb->pdebug->Leave();
}
