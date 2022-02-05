//! \file lmars.cpp
//  \brief LMARS based on Chen's 2013 paper

// C/C++ headers
#include <cstring>

// Athena++ headers
#include "../../hydro.hpp"
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../eos/eos.hpp"
#include "../../../thermodynamics/thermodynamics.hpp"
#include "../../../math/core.h" // sqr
#include "../../../globals.hpp"
#include "../../../debugger/debugger.hpp"

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu, 
                          int const ivx,
                          AthenaArray<Real> &wl, AthenaArray<Real> &wr,
                          AthenaArray<Real> &flx, const AthenaArray<Real> &dxw)
{
  Debugger *pdbg = pmy_block->pdebug;
  char str[80];
  sprintf(str, "Hydro::RiemannSolver-X%d, (k=%d,j=%d)", ivx-IVX+1, k, j);
  pdbg->Call(str);
  std::stringstream &msg = pdbg->msg;

  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Thermodynamics *pthermo = pmy_block->pthermo;
  Real rhobar, pbar, cbar, ubar, hl, hr;
  Real gamma = pmy_block->peos->GetGamma();
  Real wli[NHYDRO], wri[NHYDRO];

  for (int i = il; i <= iu; ++i) {
    // copy local variables
    for (int n = 0; n < NHYDRO; ++n) {
      wli[n] = wl(n,i);
      wri[n] = wr(n,i);
    }
    // correction for gamma
    // left
    Real fsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += wli[n]*(pthermo->GetCvRatio(n) - 1.);
      feps += wli[n]*(1./pthermo->GetMassRatio(n) - 1.);
    }
    Real kappal = 1./(gamma - 1.)*fsig/feps;

    // right
    fsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += wri[n]*(pthermo->GetCvRatio(n) - 1.);
      feps += wri[n]*(1./pthermo->GetMassRatio(n) - 1.);
    }
    Real kappar = 1./(gamma - 1.)*fsig/feps;

#if (DEBUG_LEVEL > 0)
    if (wli[IPR] + wri[IPR] < 0.) {
      msg << "### FATAL ERROR in Riemann Solver LMARS"
          << std::endl << "Pressure is negative: " << wli[IPR] << "," << wri[IPR]
          << std::endl << "At position ("
          << k << "," << j << "," << i << ") in rank " << Globals::my_rank;
      ATHENA_ERROR(msg);
    }
    if (rhobar < 0.) {
      msg << "### FATAL ERROR in Riemann Solver LMARS"
          << std::endl << "Density is negative: " << wli[IDN] << "," << wri[IDN]
          << std::endl << "At position ("
          << k << "," << j << "," << i << ") in rank " << Globals::my_rank;
      ATHENA_ERROR(msg);
    }
    if (!(kappal > 0.) || !(kappar > 0.)) {
      msg << "### FATAL ERROR in Riemann Solver LMARS"
          << std::endl << "Kappa is negative: " << kappal << "," << kappar
          << std::endl << "At position ("
          << k << "," << j << "," << i << ") in rank " << Globals::my_rank;
      ATHENA_ERROR(msg);
    }
#endif

    // enthalpy
    hl = wli[IPR]/wli[IDN]*(kappal + 1.) 
      + 0.5*(sqr(wli[ivx]) + sqr(wli[ivy]) + sqr(wli[ivz]));
    hr = wri[IPR]/wri[IDN]*(kappar + 1.)
      + 0.5*(sqr(wri[ivx]) + sqr(wri[ivy]) + sqr(wri[ivz]));

    rhobar = 0.5*(wli[IDN] + wri[IDN]);
    cbar = sqrt(0.5*(1. + (1./kappar + 1./kappal)/2.)*(wli[IPR]+ wri[IPR])/rhobar);
    pbar = 0.5*(wli[IPR] + wri[IPR]) + 0.5*(rhobar*cbar)*(wli[ivx] - wri[ivx]);
    ubar = 0.5*(wli[ivx] + wri[ivx]) + 0.5/(rhobar*cbar)*(wli[IPR] - wri[IPR]);

    Real rd = 1.;
    if (ubar > 0.) {
      // volume mixing ratio to mass mixing ratio
      for (int n = 1; n <= NVAPOR; ++n)
        rd -= wli[n];

      flx(IDN,k,j,i) = ubar*wli[IDN]*rd;
      for (int n = 1; n <= NVAPOR; ++n)
        flx(n,k,j,i) = ubar*wli[IDN]*wli[n];
      flx(ivx,k,j,i) = ubar*wli[IDN]*wli[ivx] + pbar;
      flx(ivy,k,j,i) = ubar*wli[IDN]*wli[ivy];
      flx(ivz,k,j,i) = ubar*wli[IDN]*wli[ivz];
      flx(IEN,k,j,i) = ubar*wli[IDN]*hl;
    } else {
      // volume mixing ratio to mass mixing ratio
      for (int n = 1; n <= NVAPOR; ++n)
        rd -= wri[n];

      flx(IDN,k,j,i) = ubar*wri[IDN]*rd;
      for (int n = 1; n <= NVAPOR; ++n)
        flx(n,k,j,i) = ubar*wri[IDN]*wri[n];
      flx(ivx,k,j,i) = ubar*wri[IDN]*wri[ivx] + pbar;
      flx(ivy,k,j,i) = ubar*wri[IDN]*wri[ivy];
      flx(ivz,k,j,i) = ubar*wri[IDN]*wri[ivz];
      flx(IEN,k,j,i) = ubar*wri[IDN]*hr;
    }
  }
  pdbg->Leave();
}
