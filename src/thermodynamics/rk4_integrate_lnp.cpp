#include <iostream>
#include <algorithm>
#include "thermodynamics.hpp"

void rk4_integrate_lnp(Real prim[], int isat[], Real const rcp[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma, Real dlnp)
{
  Real step[] = {0.5, 0.5, 1.};
  Real temp = prim[IDN];
  Real pres = prim[IPR];
  Real chi[4];

  for (int rk = 0; rk < 4; ++rk) {
    // reset vapor
    for (int n = 1; n < 1 + NVAPOR; ++n) {
      prim[n] += prim[n+NVAPOR];
      prim[n+NVAPOR] = 0.;
      if (isat[n] >= 0) isat[n] = 0;
    }

    for (int n = 1; n < 1 + NVAPOR; ++n) {
      int nc = n + NVAPOR;
      Real rate = VaporCloudEquilibrium(prim, n, nc, t3[nc], p3[nc], 0., beta[nc], delta[nc]);
      prim[n] -= rate;
      prim[nc] += rate;
      if ((rate > 0.) && (isat[n] >= 0)) isat[n] = 1;
    }

    // calculate tendency
    chi[rk] = dlnTdlnP(prim, isat, rcp, beta, delta, t3, gamma);
    
    // integrate over dlnp
    if (rk < 3)
      prim[IDN] = temp*exp(chi[rk]*dlnp*step[rk]);
    else
      prim[IDN] = temp*exp(1./6.*(chi[0] + 2.*chi[1] + 2.*chi[2] + chi[3])*dlnp);
    prim[IPR] = pres*exp(dlnp);
  }

  // recondensation
  for (int n = 1; n < 1 + NVAPOR; ++n) {
    int nc = n + NVAPOR;
    Real rate = VaporCloudEquilibrium(prim, n, nc, t3[n], p3[nc], 0., beta[nc], delta[nc]);
    prim[n] -= rate;
    prim[nc] += rate;
    if ((rate > 0.) && (isat[n] >= 0)) isat[n] = 1;
  }
}

void rk4_integrate_lnp_adaptive(Real prim[], int isat[], Real const rcp[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma, Real dlnp, Real ftol)
{
  Real prim1[NHYDRO], prim2[NHYDRO];
  int  isat1[1+NVAPOR], isat2[1+NVAPOR];

  for (int n = 0; n < NHYDRO; ++n) {
    prim1[n] = prim[n];
    prim2[n] = prim[n];
  }
  for (int n = 0; n < 1+ NVAPOR; ++n) {
    isat1[n] = isat[n];
    isat2[n] = isat[n];
  }

  // trail step
  rk4_integrate_lnp(prim1, isat1, rcp, beta, delta, t3, p3, gamma, dlnp);

  // refined step
  rk4_integrate_lnp(prim2, isat2, rcp, beta, delta, t3, p3, gamma, dlnp/2.);
  rk4_integrate_lnp(prim2, isat2, rcp, beta, delta, t3, p3, gamma, dlnp/2.);

  if (fabs(prim2[IDN] - prim1[IDN]) > ftol) {
    rk4_integrate_lnp_adaptive(prim, isat, rcp, beta, delta, t3, p3, gamma, dlnp/2., ftol/2.);
    rk4_integrate_lnp_adaptive(prim, isat, rcp, beta, delta, t3, p3, gamma, dlnp/2., ftol/2.);
  } else {
    for (int n = 0; n < NHYDRO; ++n)
      prim[n] = prim2[n];
    for (int n = 0; n < 1+ NVAPOR; ++n)
      isat[n] = isat2[n];
  }
}
