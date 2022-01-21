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
#include "../math/root.h"
#include "../debugger/debugger.hpp"
#include "gaussian_process.hpp"

struct SolverData {
  Thermodynamics *pthermo;
  Real **w2;
  Real dlnp;
};

Real solve_thetav(Real rdlnTdlnP, void *aux) {
  // grav parameter is not used in hydrostatic formulation, set to zero
  SolverData *pdata = static_cast<SolverData*>(aux);
  Real **w2 = pdata->w2;
  Thermodynamics *pthermo = pdata->pthermo;
  pthermo->ConstructAtmosphere(w2, pthermo->GetTemp(w2[0]), w2[0][IPR], 0., pdata->dlnp, 2, Adiabat::dry, rdlnTdlnP);
  Real p0 = 1.E5;
  Real thetav0 = PotentialTemp(w2[0], p0, pthermo)*pthermo->RovRd(w2[0]);
  Real thetav1 = PotentialTemp(w2[1], p0, pthermo)*pthermo->RovRd(w2[1]);
  return thetav1 - thetav0;
}

void update_atm_profiles(MeshBlock *pmb, int k,
    Real const *PrSample, Real const *TpSample, Real const *XpSample, int nsample, 
		std::vector<int> const& ix, Real Tstd, Real Tlen, Real Xstd, Real Xlen, Real chi)
{
  //ATHENA_LOG("update_atm_profiles");
  std::stringstream msg;
  msg << "- updating atmospheric profiles ..." << std::endl;
  Thermodynamics *pthermo = pmb->pthermo;
  Coordinates *pcoord = pmb->pcoord;
  Hydro *phydro = pmb->phydro;
  int is = pmb->is, js = pmb->js, ie = pmb->ie, je = pmb->je;

  int nlayer = ie - is + 1;
  Real *zlev = new Real [nsample];
  Real P0 = phydro->reference_pressure;
  Real H0 = phydro->scale_height;

  msg << "- sample levels: ";
  for (int i = 0; i < nsample; ++i) {
    zlev[i] = -H0*log(PrSample[i]/P0);
    msg << zlev[i] << " ";
  }
  msg << std::endl;

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
        phydro->w(n,k,j,i) = phydro->w(n,k,js,i);
	int j1 = js+1, j2 = js+2;
  Real Rd = pthermo->GetRd();

  // save perturbed T profile to model 1
	if (std::find(ix.begin(), ix.end(), 0) != ix.end()) {
		msg << "- update temperature" << std::endl;
		for (int i = is; i <= ie; ++i) {
      if (pcoord->x1v(i) < zlev[0] || pcoord->x1v(i) > zlev[nsample-1])
        continue;
			Real temp = pthermo->GetTemp(phydro->w.at(k,j1,i));
			if (temp + Tp[i-is] < 0.) Tp[i-is] = 1. - temp; // min 1K temperature
			phydro->w(IDN,k,j1,i) = phydro->w(IPR,k,j1,i)/(Rd*(temp + Tp[i-is])*
					pthermo->RovRd(phydro->w.at(k,j1,i)));
		}
	}

  // save perturbed X profile to model 2
  msg << "- update composition" << std::endl;
  for (int i = is; i <= ie; ++i) {
    Real temp = pthermo->GetTemp(phydro->w.at(k,j2,i));
    if (pcoord->x1v(i) < zlev[0] || pcoord->x1v(i) > zlev[nsample-1])
      continue;
		int ic = 0;
		for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m) {
			if (*m != 0) {
				phydro->w(*m,k,j2,i) += Xp[ic*nsample + i-is];
				phydro->w(*m,k,j2,i) = std::max(phydro->w(*m,k,j2,i), 0.);
				phydro->w(IDN,k,j2,i) = phydro->w(IPR,k,j2,i)/
					(Rd*temp*pthermo->RovRd(phydro->w.at(k,j2,i)));
				ic++;
			}
		}
  }

  // save convectively adjusted profile to model 3 (j = je)
  Real **w2, dw[1+NVAPOR];
  NewCArray(w2, 2, NHYDRO+2*NVAPOR);
  msg << "- doing convective adjustment" << std::endl;
  // save convectively adjusted profile to model 3 (j = je)
	for (int i = is+1; i <= ie; ++i) {
    if (pcoord->x1v(i) < zlev[0]) continue;
    // copy unadjusted temperature and composition profile to je
    Real temp = pthermo->GetTemp(phydro->w.at(k,j1,i));
		for (std::vector<int>::const_iterator m = ix.begin(); m != ix.end(); ++m)
			if (*m != 0) phydro->w(*m,k,je,i) = phydro->w(*m,k,j2,i);
    phydro->w(IDN,k,je,i) = phydro->w(IPR,k,je,i)/
      (Rd*temp*pthermo->RovRd(phydro->w.at(k,je,i)));

    // constant virtual potential temperature move
    for (int n = 0; n < NHYDRO; ++n)
      w2[0][n] = phydro->w(n,k,je,i-1);

    SolverData solver_data;
    solver_data.w2 = w2;
    solver_data.pthermo = pthermo;
    solver_data.dlnp = log(phydro->w(IPR,k,je,i)/phydro->w(IPR,k,je,i-1));

    Real rdlnTdlnP = 1.;
    // TODO (cli) finish convective adjustment
    /*std::cout << solve_thetav(1., &solver_data) << std::endl;
    int err = root(0.5, 2., 1.E-4, &rdlnTdlnP, solve_thetav, &solver_data);
    if (err) {
      msg << "### Root doesn't converge" << std::endl
          << solve_thetav(1., &solver_data) << " " << solve_thetav(2., &solver_data);
      ATHENA_ERROR(msg);
    }
    //msg << "- rdlnTdlnP = " << rdlnTdlnP << std::endl;*/

    pthermo->ConstructAtmosphere(w2, pthermo->GetTemp(w2[0]), w2[0][IPR], 0., solver_data.dlnp, 
        2, Adiabat::dry, rdlnTdlnP);

    // stability
    phydro->w(IDN,k,je,i) = std::min(w2[1][IDN], phydro->w(IDN,k,je,i));

    // saturation
    pthermo->SaturationSurplus(dw, phydro->w.at(k,je,i), VariableType::prim);
    for (int n = 1; n <= NVAPOR; ++n)
      if (dw[n] > 0.) phydro->w(n,k,je,i) -= dw[n];
	}
  pmb->pdebug->WriteMessage(msg.str());

  FreeCArray(w2);
  delete[] zlev;
  delete[] stdAll;
  delete[] stdSample;
  delete[] Tp;
  delete[] Xp;
}
