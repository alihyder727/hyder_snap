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
#include "../math/interpolation.h"

void calculate_fit_target(MeshBlock *pmb, Real *val, int nvalue,
    int k, int j, bool differential)
{
  //ATHENA_LOG("calculate_fit_target");
  std::stringstream &msg = pmb->pdebug->msg;
  std::vector<Direction> out_dir = pmb->prad->GetOutgoingRays();
  int ndir = out_dir.size();

  Real *mus = new Real [ndir];
  Real *tbs = new Real [ndir];

  pmb->pdebug->Call("calculate_fit_target");
  msg << "- model " << j - pmb->js << std::endl;
  for (int i = 0; i < ndir; ++i)
    mus[i] = out_dir[i].mu;

  // 11. log likelihood
  RadiationBand *pband = pmb->prad->pband;
  int i = 0;
  while (pband != NULL) {
    // brightness temperatures
    val[i*2] = pband->btoa(0,k,j);

    // limb darkening
    for (int n = 0; n < ndir; ++n)
      tbs[n] = pband->btoa(n,k,j);
    
    Real tb45 = interp1(cos(45./180.*M_PI), tbs, mus, ndir);
    val[i*2+1] = (tbs[0] - tb45)/tbs[0]*100.;

    if (differential) {
      // brightness temperatures
      val[i*2] -= pband->btoa(0,k,pmb->js);

      // limb darkening
      for (int n = 0; n < ndir; ++n)
        tbs[n] = pband->btoa(n,k,pmb->js);

      tb45 = interp1(cos(45./180.*M_PI), tbs, mus, ndir);
      val[i*2+1] -= (tbs[0] - tb45)/tbs[0]*100.;
    }

    i++;
    pband = pband->next;
  }

  if (nvalue != 2*i) {
    msg << "### FATAL ERROR in calculate_fit_target" << std::endl
        << "Lengths of forward vector and target vector do not match";
    ATHENA_ERROR(msg);
  }

  msg << "- foward model results: ";
  for (int i = 0; i < nvalue; ++i)
    msg << std::setprecision(5) << val[i] << " ";
  msg << std::endl;

  delete[] mus;
  delete[] tbs;

  pmb->pdebug->Leave();
}
