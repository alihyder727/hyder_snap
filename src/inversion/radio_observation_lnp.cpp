// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../radiation/radiation.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/root.h"
#include "../math/linalg.h"
#include "inversion.hpp"
#include "gaussian_process.hpp"
#include "radio_data.hpp"

struct RootData {
  Real rho0;
  Real Tv0;
  Real dz;
  Real Rd;
  Real grav;
  Real chi;
  Real p0;
};

Real junomwr_root_func(Real Tv2, void *aux) {
  RootData *pd = static_cast<RootData*>(aux);
  Real& rho0 = pd->rho0;
  Real& Tv0 = pd->Tv0;
  Real& dz = pd->dz;
  Real& Rd = pd->Rd;
  Real& grav = pd->grav;
  Real& chi = pd->chi;
  Real& p0 = pd->p0;

  Real rho2 = rho0*(Tv0 - grav*dz/(2.*Rd))/(Tv2 + grav*dz/(2.*Rd));
  Real p2 = rho2*Rd*Tv2;
  if (rho2 > 0.)
    return Tv2 - Tv0*pow(p2/p0, chi);
  else
    return (Tv2 + Tv0)*log(Tv2/Tv0) + 2.*chi*grav*dz/Rd;
}

Real RadioObservationLnProb(Real *par, Real *val, int ndim, int nvalue, void *obj)
{
  ATHENA_LOG("RadioObservationLnProb");
  // 1. load fitting data 
  RadioData *pobj = static_cast<RadioData*>(obj);

  // 2. forward modeling
  MeshBlock *pmb = pobj->pmy_block;
  Hydro *phydro = pmb->phydro;
  Radiation *prad = pmb->prad;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  int nsample = ndim/3;
  Real Z0 = 0., P0 = 1.E5;

  Real *PrSample = new Real [nsample];
  Real *TpSample = new Real [nsample];
  Real *XpSample = new Real [nsample];

  // 0..nsample-1 : sample pressure , the first one is fixed
  // nsample..2*nsample-1 : sample temperature
  // 2*nsample..3*nsample-1 : sample composition
  std::cout << "- Parameters: ";
  for (int i = 0; i < nsample; ++i) {
    PrSample[i] = par[i]*1.E5; // bar -> pa
    TpSample[i] = par[nsample+i];
    XpSample[i] = par[2*nsample+i]/1.E3;  // g/kg -> kg/kg
    std::cout << par[i] << " " << par[nsample+i] << " " << par[2*nsample+i] << " ";
  }
  std::cout << std::endl;

  update_atm_profiles(pmb, PrSample, TpSample, XpSample, nsample, 
      pobj->ix, pobj->Tstd, pobj->Tlen, pobj->Xstd, pobj->Xlen, P0, Z0);

  // calculate radiation for updated profiles
  for (int j = js; j <= je; ++j)
    prad->CalculateRadiances(phydro->w, 0., ks, j, is, ie+1);

  // 3. log probability
  calculate_fit_target(pmb, js+1, val, nvalue);
  Eigen::VectorXd misfit(nvalue);

  for (int m = 0; m < nvalue; ++m)
    misfit(m) = val[m] - pobj->target(m);

  Real lnprob = -0.5*misfit.transpose()*pobj->icov*misfit;

  // prior probability
  Real *zlev = new Real [nsample];
  Real *stdSample = new Real [nsample];
  for (int i = 0; i < nsample; ++i)
    zlev[i] = Z0 - phydro->scale_height*log(PrSample[i]/P0);

  for (int i = 0; i < nsample; ++i)
    stdSample[i] = pobj->Xstd;
  Real lnp1 = gp_lnprior(SquaredExponential, XpSample, zlev, stdSample, nsample, pobj->Xlen);

  for (int i = 0; i < nsample; ++i)
    stdSample[i] = pobj->Tstd;
  Real lnp2 = gp_lnprior(SquaredExponential, TpSample, zlev, stdSample, nsample, pobj->Tlen);

  std::cout << "- Log probabilities: ";
  std::cout << lnp1 << " " << lnp2 << " " << lnprob << std::endl;

  lnprob += lnp1 + lnp2;

  delete[] PrSample;
  delete[] TpSample;
  delete[] XpSample;
  delete[] zlev;
  delete[] stdSample;

  return lnprob;
}
