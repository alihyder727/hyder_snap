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
#include "../math/linalg.h"

void calculate_fit_target(MeshBlock *pmb, int j, Real *val, int nvalue)
{
  ATHENA_LOG("calculate_fit_target");
  std::stringstream msg;
  std::vector<Direction> out_dir = pmb->prad->GetOutgoingRays();
  int nangle = 0;
  std::cout << "* Model id: " << j - pmb->js << std::endl;
  std::cout << "* Emission angles used: ";
  for (int i = 0; i < out_dir.size(); ++i)
    if (out_dir[i].mu >= cos(45./180.*M_PI)) {
      std::cout << acos(out_dir[i].mu)/M_PI*180 << " ";
      nangle++;
    }
  std::cout << std::endl;

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
      b[n] = pband->btoa(n,pmb->ks,j);
    leastsq(B, b, nangle, 3);

    for (int n = 0; n < 3; ++n)
      val[i*3+n] = b[n];

    i++;
    pband = pband->next;
  }

  if (nvalue != 3*i) {
    msg << "### FATAL ERROR in calculate_fit_target" << std::endl
        << "Lengths of forward vector and target vector do not match";
    ATHENA_ERROR(msg);
  }

  std::cout << "* Foward model results: ";
  for (int i = 0; i < nvalue; ++i)
    std::cout << std::setprecision(4) << val[i] << " ";
  std::cout << std::endl;

  FreeCArray(B);
  delete[] b;
}

