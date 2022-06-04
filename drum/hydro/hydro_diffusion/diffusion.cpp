//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "hydro_diffusion.hpp"

//---------------------------------------------------------------------------------------
// Calculate isotropic diffusion
// diffusion based on concentration gradient
//! \todo may not work for non-Cartesian geometry

void HydroDiffusion::DiffusionFluxIso(
    const AthenaArray<Real> &prim,
    const AthenaArray<Real> &cons, AthenaArray<Real> *dfuflx) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  AthenaArray<Real> &x1flux = dfuflx[X1DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real nuf, denf, dqdx, dqdy, dqdz;

  // i-direction
  jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) {
      if (!f3) // 2D
        jl = js - 1, ju = je + 1, kl = ks, ku = ke;
      else // 3D
        jl = js - 1, ju = je + 1, kl = ks - 1, ku = ke + 1;
    }
  }
  // x1flux = nu*rho*dqdx/schmidt
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        nuf = 0.5*(nu(DiffProcess::iso,k,j,i) + nu(DiffProcess::iso,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
        for (int n = 1; n <= NVAPOR; ++n) {
          dqdx = (prim(n,k,j,i) - prim(n,k,j,i-1))/pco_->dx1v(i-1);
          x1flux(n,k,j,i) -= nuf*denf*dqdx/schmidt_;
        }
      }
    }
  }

  // j-direction
  il = is, iu = ie, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (!f3) // 2D
      il = is - 1, iu = ie + 1, kl = ks, ku = ke;
    else // 3D
      il = is - 1, iu = ie + 1, kl = ks - 1, ku = ke + 1;
  }
  if (f2) { // 2D or 3D
    AthenaArray<Real> &x2flux = dfuflx[X2DIR];
    for (int k=kl; k<=ku; ++k) {
      for (int j=js; j<=je+1; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          nuf = 0.5*(nu(DiffProcess::iso,k,j,i) + nu(DiffProcess::iso,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j-1,i));
          for (int n = 1; n <= NVAPOR; ++n) {
            dqdy = (prim(n,k,j,i) - prim(n,k,j-1,i))/
              pco_->h2v(i)/pco_->h2v(i)/pco_->dx2v(j-1);
            x2flux(n,k,j,i) -= nuf*denf*dqdy/schmidt_;
          }
        }
      }
    }
  } // zero flux for 1D

  // k-direction
  il = is, iu = ie, jl = js, ju = je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) // 2D or 3D
      il = is - 1, iu = ie + 1, jl = js - 1, ju = je + 1;
    else // 1D
      il = is - 1, iu = ie + 1;
  }
  if (f3) { // 3D
    AthenaArray<Real> &x3flux = cndflx[X3DIR];
    for (int k=ks; k<=ke+1; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          nuf = 0.5*(nu(DiffProcess::iso,k,j,i) + nu(DiffProcess::iso,k-1,j,i));
          denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k-1,j,i));
          for (int n = 1; n <= NVAPOR; ++n) {
            dqdz = (prim(n,k,j,i) - prim(n,k-1,j,i))/
              pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
            x3flux(n,k,j,i) -= nuf*denf*dqdz/schmidt_;
          }
        }
      }
    }
  } // zero flux for 1D/2D
  return;
}


//---------------------------------------------------------------------------------------
// Calculate anisotropic diffusion

void HydroDiffusion::DiffusionFluxAniso(
    const AthenaArray<Real> &p,
    const AthenaArray<Real> &c, AthenaArray<Real> *flx) {
  return;
}
