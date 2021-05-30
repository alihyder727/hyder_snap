/** @file constuct_atmosphere.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Wednesday May 26, 2021 18:45:55 UTC
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>
#include <cstdlib>

// Athena++ headers
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "thermodynamics.hpp"
#include "thermodynamic_funcs.hpp"

void Thermodynamics::ConstructAtmosphere(Real **w, Real Ts, Real Ps,
  Real grav, Real dz, int len, Adiabat method, Real dTdz) const
{
  Real gamma = pmy_block->peos->GetGamma();

  // hydro + liquid + ice
  Real q1[NHYDRO+2*NVAPOR];
  PrimitiveToChemical(q1, w[0]);
  for (int n = NHYDRO; n < NHYDRO+2*NVAPOR; ++n) q1[n] = 0.;
  q1[IDN] = Ts;
  q1[IPR] = Ps;

  Real rcp[1+3*NVAPOR];  // molar cp ratio
  for (int n = 0; n < 1+3*NVAPOR;  ++n)
    rcp[n] = cp_ratios_[n]*mu_ratios_[n];

  int isat[1+NVAPOR];
  for (int n = 0; n <= NVAPOR; ++n)
    isat[n] = 0;

  // equilibrate with liquid (iv+NVAPOR) or ice (iv+2*NVAPOR)
  Real xg = 1.; // total moles in gas phase
  for (int iv = 1; iv <= NVAPOR; ++iv) {
    int nc = q1[IDN] > t3_[iv] ? iv + NVAPOR : iv + 2*NVAPOR;
    int ic = NHYDRO - NVAPOR + nc - 1;
    Real rate = VaporCloudEquilibrium(q1, iv, ic, t3_[iv], p3_[iv], 
        0., beta_[nc], delta_[nc]);
    q1[iv] -= rate;
    q1[ic] += rate;
    if (rate > 0.) isat[iv] = 1;
    xg -= rate;
  }

  // re-calculate vapor molar mixing ratios
  for (int iv = 1; iv <= NVAPOR; ++iv) q1[iv] /= xg;
  // vapor
  ChemicalToPrimitive(w[0], q1);
  // set back
  for (int iv = 1; iv <= NVAPOR; ++iv) q1[iv] *= xg;
  // clouds, density
  // $\frac{\rho_{ij}}{\rho_1} = \frac{\hat{q}_{ij}}{\hat{q}_1}\frac{\mu_i}{\mu_1}
  for (int n = 0; n < 2*NVAPOR; ++n)
    w[0][NHYDRO+n] = (w[0][0]*w[0][1])*(q1[NHYDRO+n]/q1[1])
                    *(mu_ratios_[1+NVAPOR+n]/mu_ratios_[1]);

  for (int i = 1; i < len; ++i) {
    // RK4 integration 
    rk4_integrate_z_adaptive(q1, isat, rcp, mu_ratios_, beta_, delta_, t3_, p3_, gamma,
      grav/Rd_, dz, ftol_, (int)method, dTdz);

    xg = 1.;
    for (int n = 0; n < 2*NVAPOR; ++n) xg -= q1[NHYDRO+n];
    for (int iv = 1; iv <= NVAPOR; ++iv) q1[iv] /= xg;
    // vapor
    ChemicalToPrimitive(w[i], q1);
    // set back
    for (int iv = 1; iv <= NVAPOR; ++iv) q1[iv] *= xg;
    // clouds
    for (int n = 0; n < 2*NVAPOR; ++n)
      w[i][NHYDRO+n] = (w[i][0]*w[i][1])*(q1[NHYDRO+n]/q1[1])
                      *(mu_ratios_[1+NVAPOR+n]/mu_ratios_[1]);
  }
}
