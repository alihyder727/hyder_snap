/** @file update_atm_profiles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Nov 18, 2021 20:28:43 EST
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../radiation/radiation.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "gaussian_process.hpp"

void update_atm_profiles(MeshBlock *pmb,
    Real *PrSample, Real *TpSample, Real *XpSample, int nsample, int ix,
    Real Tstd, Real Tlen, Real Xstd, Real Xlen, Real P0, Real Z0 = 0.)
{
  ATHENA_LOG("update_atm_profiles");
  Thermodynamics *pthermo = pmb->pthermo;
  Coordinates *pcoord = pmb->pcoord;
  Hydro *phydro = pmb->phydro;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  int nlayer = ie - is + 1;
  Real *zlev = new Real [nsample];

  std::cout << "* Sample levels: ";
  for (int i = 0; i < nsample; ++i) {
    zlev[i] = Z0 - phydro->scale_height*log(PrSample[i]/P0);
    std::cout << zlev[i] << " ";
  }
  std::cout << std::endl;

  // calculate the covariance matrix of T
  Real *stdAll = new Real [nlayer];
  Real *stdSample = new Real [nsample];
  Real *Tp = new Real [nlayer];
  Real *Xp = new Real [nlayer];

  // calculate perturbed T profile
  for (int i = is; i <= ie; ++i)
    stdAll[i-is] = Tstd;
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Tstd;

  gp_predict(SquaredExponential, Tp, &pcoord->x1v(is), stdAll, nlayer,
    TpSample, zlev, stdSample, nsample, Tlen);

  // fix boundary condition
  int ib = 0, it = nlayer - 1;
  while ((Tp[ib+1] < Tp[ib]) && (ib < nlayer)) ib++;
  for (int i = 0; i < ib; ++i) Tp[i] = Tp[ib]; 
  while ((Tp[it-1] > Tp[it]) && (it >= 0)) it--;
  for (int i = it+1; i < nlayer; ++i) Tp[i] = Tp[it]; 

  // calculate perturbed X profile
  for (int i = is; i <= ie; ++i)
    stdAll[i-is] = Xstd;
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Xstd;
    
  gp_predict(SquaredExponential, Xp, &pcoord->x1v(is), stdAll, nlayer,
    XpSample, zlev, stdSample, nsample, Xlen);

  // copy baseline js -> js+1 .. je
  for (int n = 0; n < NHYDRO; ++n)
    for (int j = js+1; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        phydro->w(n,j,i) = phydro->w(n,js,i);

  // save perturbed X profile
  std::cout << "* Calculate Tb if only X was perturbed" << std::endl;
  Real rho, Rd = pthermo->GetRd();
  int j1 = js+1, j2 = js+2, j3 = js+3;
  for (int i = is; i <= ie; ++i) {
    Real temp = pthermo->GetTemp(phydro->w.at(j1,i));
    phydro->w(ix,j1,i) += Xp[i-is];
    phydro->w(ix,j1,i) = std::max(phydro->w(ix,j1,i), 0.);
    phydro->w(IDN,j1,i) = phydro->w(IPR,j1,i)/
      (Rd*temp*pthermo->RovRd(phydro->w.at(j1,i)));
  }

  // save perturbed T profile
  std::cout << "* Calculate Tb if only T was perturbed" << std::endl;
  for (int i = is; i <= ie; ++i) {
    Real temp = pthermo->GetTemp(phydro->w.at(j2,i));
    if (temp + Tp[i-is] < 0.) Tp[i-is] = 1. - temp; // min 1K temperature
    phydro->w(IDN,j2,i) = phydro->w(IPR,j2,i)/(Rd*(temp + Tp[i-is])*
        pthermo->RovRd(phydro->w.at(j2,i)));
  }

  // convective adjustment

  delete[] zlev;
  delete[] stdAll;
  delete[] stdSample;
  delete[] Tp;
  delete[] Xp;
}
