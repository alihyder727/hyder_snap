/** @file update_atm_profiles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Thursday Nov 18, 2021 20:28:43 EST
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>
#include <algorithm>

// Athena++ headers
#include "../athena.hpp"
#include "../mesh/mesh.hpp"
#include "../radiation/radiation.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../utils/utils.hpp"
#include "../thermodynamics/thermodynamic_funcs.hpp"
#include "gaussian_process.hpp"

void update_atm_profiles(MeshBlock *pmb,
    Real const *PrSample, Real const *TpSample, Real const *XpSample, int nsample, 
		std::vector<int> const& ix, Real Tstd, Real Tlen, Real Xstd, Real Xlen, Real chi)
{
  ATHENA_LOG("update_atm_profiles");
  Thermodynamics *pthermo = pmb->pthermo;
  Coordinates *pcoord = pmb->pcoord;
  Hydro *phydro = pmb->phydro;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  int nlayer = ie - is + 1;
  Real *zlev = new Real [nsample];
  Real P0 = phydro->reference_pressure;
  Real H0 = phydro->scale_height;

  std::cout << "* Sample levels: ";
  for (int i = 0; i < nsample; ++i) {
    zlev[i] = -H0*log(PrSample[i]/P0);
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
    stdAll[i-is] = Tstd*pow(exp(pcoord->x1v(i)/H0), chi);
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Tstd*pow(exp(zlev[i]/H0), chi);

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
    stdAll[i-is] = Xstd*pow(exp(pcoord->x1v(i)/H0), chi);
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Xstd*pow(exp(zlev[i]/H0), chi);
    
  gp_predict(SquaredExponential, Xp, &pcoord->x1v(is), stdAll, nlayer,
    XpSample, zlev, stdSample, nsample, Xlen);

  // copy baseline js -> js+1 .. je
  for (int n = 0; n < NHYDRO; ++n)
    for (int j = js+1; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        phydro->w(n,j,i) = phydro->w(n,js,i);
	int j1 = js+1, j2 = js+2;
  Real Rd = pthermo->GetRd();

  // save perturbed T profile to model 1
	if (std::find(ix.begin(), ix.end(), 0) != ix.end()) {
		std::cout << "* Update temperature" << std::endl;
		for (int i = is; i <= ie; ++i) {
      if (pcoord->x1v(i) < zlev[0] || pcoord->x1v(i) > zlev[nsample-1])
        continue;
			Real temp = pthermo->GetTemp(phydro->w.at(j1,i));
			if (temp + Tp[i-is] < 0.) Tp[i-is] = 1. - temp; // min 1K temperature
			phydro->w(IDN,j1,i) = phydro->w(IPR,j1,i)/(Rd*(temp + Tp[i-is])*
					pthermo->RovRd(phydro->w.at(j1,i)));
		}
	}

  // save perturbed X profile to model 2
  std::cout << "* Update composition" << std::endl;
  for (int i = is; i <= ie; ++i) {
    Real temp = pthermo->GetTemp(phydro->w.at(j2,i));
    if (pcoord->x1v(i) < zlev[0] || pcoord->x1v(i) > zlev[nsample-1])
      continue;
		int ic = 0;
		for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m) {
			if (*m != 0) {
				phydro->w(*m,j2,i) += Xp[ic*nsample + i-is];
				phydro->w(*m,j2,i) = std::max(phydro->w(*m,j2,i), 0.);
				phydro->w(IDN,j2,i) = phydro->w(IPR,j2,i)/
					(Rd*temp*pthermo->RovRd(phydro->w.at(j2,i)));
				ic++;
			}
		}
  }

  // save convectively adjusted profile to model 3 (j = je)
  Real **w2, dw[1+NVAPOR];
  NewCArray(w2, 2, NHYDRO+2*NVAPOR);
  std::cout << "* Convective adjustment" << std::endl;
  // save convectively adjusted profile to model 3 (j = je)
	for (int i = is+1; i <= ie; ++i) {
    if (pcoord->x1v(i) < zlev[0]) continue;
    // copy unadjusted temperature and composition profile to je
    Real temp = pthermo->GetTemp(phydro->w.at(j1,i));
		for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m)
			if (*m != 0) phydro->w(*m,je,i) = phydro->w(*m,j2,i);
    phydro->w(IDN,je,i) = phydro->w(IPR,je,i)/
      (Rd*temp*pthermo->RovRd(phydro->w.at(je,i)));

    // adiabatic move
    for (int n = 0; n < NHYDRO; ++n)
      w2[0][n] = phydro->w(n,je,i-1);
    temp = pthermo->GetTemp(w2[0]);
    Real pres = w2[0][IPR];
    Real dlnp = log(phydro->w(IPR,je,i)/pres);
    // grav parameter is not used in hydrostatic formulation, set to zero
    pthermo->ConstructAtmosphere(w2, temp, pres, 0., dlnp, 2, Adiabat::pseudo, 1.);
  
    // stability
    phydro->w(IDN,je,i) = std::min(w2[1][IDN], phydro->w(IDN,je,i));

    // saturation
    pthermo->SaturationSurplus(dw, phydro->w.at(je,i), VariableType::prim);
    for (int n = 1; n <= NVAPOR; ++n)
      if (dw[n] > 0.) phydro->w(n,je,i) -= dw[n];
	}
  FreeCArray(w2);

  delete[] zlev;
  delete[] stdAll;
  delete[] stdSample;
  delete[] Tp;
  delete[] Xp;
}
