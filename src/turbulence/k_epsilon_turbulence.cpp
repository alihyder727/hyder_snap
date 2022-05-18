//! \file k_epsilon_turbulence.cpp
//  \brief implement K-Epsilon turbulence

// C headers

// C++ headers
#include <algorithm>   // min,max

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"   // reapply floors to face-centered reconstructed states
#include "../hydro/hydro.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "turbulence_model.hpp"

KEpsilonTurbulence::KEpsilonTurbulence(MeshBlock *pmb, ParameterInput *pin):
  TurbulenceModel(pmb, pin, 2)
{
  cmu_ = pin->GetOrAddReal("turbulence", "kepsilon.cmu", 0.09);
  c1_ = pin->GetOrAddReal("turbulence", "kepsilon.c1", 1.44);
  c2_ = pin->GetOrAddReal("turbulence", "kepsilon.c2", 1.92);
  sigk_ = pin->GetOrAddReal("turbulence", "kepsilon.sigk", 1.0);
  sige_ = pin->GetOrAddReal("turbulence", "kepsilon.sige", 1.3);
}

inline Real Laplace_(AthenaArray<Real> const& mut, AthenaArray<Real> const& v,
  int k, int j, int i, Coordinates *pcoord)
{
  Real result = 0.;
  Mesh *pm = pcoord->pmy_block->pmy_mesh;

  // first dimension
  Real mut_p = (mut(k,j,i+1) + mut(k,j,i))/2.;
  Real mut_m = (mut(k,j,i) + mut(k,j,i-1))/2.;
  Real gradv_p = (v(k,j,i+1) - v(k,j,i))/(pcoord->x1v(i+1) - pcoord->x1v(i));
  Real gradv_m = (v(k,j,i) - v(k,j,i-1))/(pcoord->x1v(i) - pcoord->x1v(i-1));
  result += (mut_p*gradv_p - mut_m*gradv_m)/pcoord->dx1f(i);

  // second dimension
  if (pm->f2) {
    mut_p = (mut(k,j+1,i) + mut(k,j,i))/2.;
    mut_m = (mut(k,j,i) + mut(k,j-1,i))/2.;
    gradv_p = (v(k,j+1,i) - v(k,j,i))/(pcoord->x2v(j+1) - pcoord->x2v(i));
    gradv_m = (v(k,j,i) - v(k,j-1,i))/(pcoord->x2v(j) - pcoord->x2v(j-1));
    result += (mut_p*gradv_p - mut_m*gradv_m)/pcoord->dx2f(j);
  }

  // third dimension
  if (pm->f3) {
    mut_p = (mut(k+1,j,i) + mut(k,j,i))/2.;
    mut_m = (mut(k,j,i) + mut(k-1,j,i))/2.;
    gradv_p = (v(k+1,j,i) - v(k,j,i))/(pcoord->x3v(k+1) - pcoord->x3v(k));
    gradv_m = (v(k,j,i) - v(k-1,j,i))/(pcoord->x3v(k) - pcoord->x3v(k-1));
    result += (mut_p*gradv_p - mut_m*gradv_m)/pcoord->dx3f(k);
  }

  return result;
}

inline Real ShearProduction_(AthenaArray<Real> const& w,
  int k, int j, int i, Coordinates *pcoord)
{
  Real result = 0.;
  Mesh *pm = pcoord->pmy_block->pmy_mesh;
  Real gradv, gradv1, gradv2, gradv3, gradv4;

  // first dimension
  gradv = (w(IM1,k,j,i+1) - w(IM1,k,j,i-1))/(pcoord->x1v(i+1) - pcoord->x1v(i-1));
  result += 2*gradv*gradv;

  // second dimension
  if (pm->f2) {
    gradv = (w(IM2,k,j+1,i) - w(IM2,k,j-1,i))/(pcoord->x2v(j+1) - pcoord->x2v(j-1));
    gradv1 = (w(IM1,k,j+1,i) - w(IM1,k,j-1,i))/(pcoord->x2v(j+1) - pcoord->x2v(j-1));
    gradv2 = (w(IM2,k,j,i+1) - w(IM2,k,j,i-1))/(pcoord->x1v(i+1) - pcoord->x1v(i-1));
    result += 2*gradv*gradv + (gradv1 + gradv2)*(gradv1 + gradv2);
  }

  // thrid dimension
  if (pm->f3) {
    gradv = (w(IM3,k+1,j,i) - w(IM3,k-1,j,i))/(pcoord->x3v(k+1) - pcoord->x3v(k-1));
    gradv1 = (w(IM1,k+1,j,i) - w(IM1,k-1,j,i))/(pcoord->x3v(k+1) - pcoord->x3v(k-1));
    gradv2 = (w(IM3,k,j,i+1) - w(IM3,k,j,i-1))/(pcoord->x1v(i+1) - pcoord->x1v(i-1));
    gradv3 = (w(IM2,k+1,j,i) - w(IM2,k-1,j,i))/(pcoord->x3v(k+1) - pcoord->x3v(k-1));
    gradv4 = (w(IM3,k,j+1,i) - w(IM3,k,j-1,i))/(pcoord->x2v(j+1) - pcoord->x2v(j-1));
    result += 2*gradv*gradv + (gradv1 + gradv2)*(gradv1 + gradv2) 
      + (gradv3 + gradv4)*(gradv3 + gradv4);
  }

  return result;
}

void KEpsilonTurbulence::driveTurbulence(AthenaArray<Real> &s,
  AthenaArray<Real> const& r, AthenaArray<Real> const& w, Real dt)
{
  MeshBlock *pmb = pmy_block;
  int ks = pmb->ks, js = pmb->js, is = pmb->is;
  int ke = pmb->ke, je = pmb->je, ie = pmb->ie;

  AthenaArray<Real> eps, tke;
  eps.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(r),4,0,1);
  tke.InitWithShallowSlice(const_cast<AthenaArray<Real>&>(r),4,1,1);

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        // shear production
        Real shear = ShearProduction_(w,k,j,i,pmb->pcoord);

        // turbulent dissipation, de/dt, eq2.2-1
        s(0,k,j,i) += (Laplace_(mut,eps,k,j,i,pmb->pcoord)/sige_
          + c1_*mut(k,j,i)*eps(k,j,i)/tke(k,j,i)*shear
          - c2_*eps(k,j,i)*eps(k,j,i)/tke(k,j,i)*w(IDN,k,j,i))*dt;

        // turbulent kinetic energy, dk/dt, eq 2.2-2
        s(1,k,j,i) += (Laplace_(mut,tke,k,j,i,pmb->pcoord)/sigk_
          + mut(k,j,i)*shear
          - eps(k,j,i)*w(IDN,k,j,i))*dt;

        // dynamic turbulent viscosity, mu_t = c_mu*k^2/epsilon, eq2.2-3
        if (eps(k,j,i) > 0.)
          mut(k,j,i) = cmu_*w(IDN,k,j,i)*tke(k,j,i)*tke(k,j,i)/eps(k,j,i);
        else
          mut(k,j,i) = 0.;
      }
}
