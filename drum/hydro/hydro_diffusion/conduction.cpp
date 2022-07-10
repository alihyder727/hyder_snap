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
#include "../../thermodynamics/thermodynamics.hpp"
#include "../../eos/eos.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_diffusion.hpp"

//---------------------------------------------------------------------------------------
// Calculate isotropic thermal conduction
// diffusion based on enthalpy gradient
// !\todo may not work for non-Cartesian geometry

void HydroDiffusion::ThermalFluxIso(
    const AthenaArray<Real> &prim,
    const AthenaArray<Real> &cons, AthenaArray<Real> *cndflx) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  AthenaArray<Real> &x1flux = cndflx[X1DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;
  Real kappaf, denf, dhdx, dhdy, dhdz, dTdx, dTdy, dTdz;
  Real gamma = pmb_->peos->GetGamma();
  Thermodynamics *pthermo = pmb_->pthermo;
  Hydro *phydro = pmb_->phydro;

  // i-direction
  il = is, iu = ie+1, jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) {
      if (!f3) // 2D
        jl = js - 1, ju = je + 1, kl = ks, ku = ke;
      else // 3D
        jl = js - 1, ju = je + 1, kl = ks - 1, ku = ke + 1;
    }
  }

  //! \todo there gonna be a better way to do this
  if (pmb_->pbval->isPhysicalBoundary(BoundaryFace::inner_x1)) il++;
  if (pmb_->pbval->isPhysicalBoundary(BoundaryFace::outer_x1)) iu--;
  Real Rd = pthermo->GetRd();

  // x1flux = kappa/gamma*rho*dhdx
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        kappaf = 0.5*(kappa(DiffProcess::iso,k,j,i) + kappa(DiffProcess::iso,k,j,i-1));
        denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j,i-1));
        //Real temp = pthermo->GetTemp(prim.at(k,j,i));
        //Real h1 = pthermo->getSpecificEnthalpy(prim.at(k,j,i));
        //Real h2 = pthermo->getSpecificEnthalpy(prim.at(k,j,i-1));
        //dhdx = (h1 - h2)/pco_->dx1v(i-1)
        Real T1 = pthermo->GetTemp(prim.at(k,j,i));
        Real T2 = pthermo->GetTemp(prim.at(k,j,i-1));
        dTdx = (T1 - T2)/pco_->dx1v(i-1);
        x1flux(k,j,i) -= kappaf*denf*dTdx*Rd/(gamma - 1.);
      }
    }
  }

  // j-direction
  il = is, iu = ie, jl = js, ju = je+1, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (!f3) // 2D
      il = is - 1, iu = ie + 1, kl = ks, ku = ke;
    else // 3D
      il = is - 1, iu = ie + 1, kl = ks - 1, ku = ke + 1;
  }
  if (pmb_->pbval->isPhysicalBoundary(BoundaryFace::inner_x2)) jl++;
  if (pmb_->pbval->isPhysicalBoundary(BoundaryFace::outer_x2)) ju--;

  if (f2) { // 2D or 3D
    AthenaArray<Real> &x2flux = cndflx[X2DIR];
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(DiffProcess::iso,k,j,i) + kappa(DiffProcess::iso,k,j-1,i));
          denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k,j-1,i));
          //Real h1 = pthermo->getSpecificEnthalpy(prim.at(k,j,i));
          //Real h2 = pthermo->getSpecificEnthalpy(prim.at(k,j-1,i));
          //dhdy = (h1 - h2)/pco_->h2v(i)/pco_->h2v(i)/pco_->dx2v(j-1) - phydro->hsrc.GetG2();
          //x2flux(k,j,i) -= kappaf*denf*dhdy/gamma;
          Real T1 = pthermo->GetTemp(prim.at(k,j,i));
          Real T2 = pthermo->GetTemp(prim.at(k,j-1,i));
          dTdy = (T1 - T2)/pco_->h2v(i)/pco_->h2v(i)/pco_->dx2v(j-1);
          x2flux(k,j,i) -= kappaf*denf*dTdy*Rd/(gamma - 1.);
        }
      }
    }
  } // zero flux for 1D

  // k-direction
  il = is, iu = ie, jl = js, ju = je, kl = ks, ku = ke+1;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) // 2D or 3D
      il = is - 1, iu = ie + 1, jl = js - 1, ju = je + 1;
    else // 1D
      il = is - 1, iu = ie + 1;
  }
  if (pmb_->pbval->isPhysicalBoundary(BoundaryFace::inner_x3)) kl++;
  if (pmb_->pbval->isPhysicalBoundary(BoundaryFace::outer_x3)) ku--;

  if (f3) { // 3D
    AthenaArray<Real> &x3flux = cndflx[X3DIR];
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          kappaf = 0.5*(kappa(DiffProcess::iso,k,j,i) + kappa(DiffProcess::iso,k-1,j,i));
          denf = 0.5*(prim(IDN,k,j,i) + prim(IDN,k-1,j,i));
          //Real h1 = pthermo->getSpecificEnthalpy(prim.at(k,j,i));
          //Real h2 = pthermo->getSpecificEnthalpy(prim.at(k-1,j,i));
          //dhdz = (h1 - h2)/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j) -
          //  phydro->hsrc.GetG3();
          //x3flux(k,j,i) -= kappaf*denf*dhdz/gamma;
          Real T1 = pthermo->GetTemp(prim.at(k,j,i));
          Real T2 = pthermo->GetTemp(prim.at(k-1,j,i));
          dTdz = (T1 - T2)/pco_->dx3v(k-1)/pco_->h31v(i)/pco_->h32v(j);
          x3flux(k,j,i) -= kappaf*denf*dTdz*Rd/(gamma - 1.);
        }
      }
    }
  } // zero flux for 1D/2D
  return;
}


//---------------------------------------------------------------------------------------
// Calculate anisotropic thermal conduction

void HydroDiffusion::ThermalFluxAniso(
    const AthenaArray<Real> &p,
    const AthenaArray<Real> &c, AthenaArray<Real> *flx) {
  return;
}


//----------------------------------------------------------------------------------------
// constant viscosity

void ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &bcc,
                     int is, int ie, int js, int je, int ks, int ke) {
  if (phdif->kappa_iso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(HydroDiffusion::DiffProcess::iso,k,j,i) = phdif->kappa_iso;
      }
    }
  }
  if (phdif->kappa_aniso > 0.0) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          phdif->kappa(HydroDiffusion::DiffProcess::aniso,k,j,i) = phdif->kappa_aniso;
      }
    }
  }
  return;
}
