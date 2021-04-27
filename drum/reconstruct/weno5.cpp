//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file weno5_simple.cpp
//  \brief  WENO5 interpolation

// Athena++ headers
#include "reconstruction.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "interpolation.hpp"
//#include "interp_weno5.hpp"
//#include "interp_weno3.hpp"


//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X1()
//  \brief 

void Reconstruction::Weno5X1(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  MeshBlock *pmb = pmy_block_;
  Real v[5], scale;
  for (int n=0; n<NMASS; ++n) {
    scale = 1.E-16;
#pragma omp simd
    for (int m = -2; m <= 2; ++m) scale += w(n,k,j,il+m);
    for (int i=il; i<=iu; ++i) {
#pragma omp simd
      for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k,j,i+m)/(5.*scale);
      wl(n,i+1) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
      wr(n,i  ) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
      scale += w(n,k,j,i+3) - w(n,k,j,i-2);
    }
  }

  int ng1 = 0, ng2 = 0;
  if (pmb->pbval->block_bcs[inner_x1] != BoundaryFlag::block)
    ng1 = NGHOST;
  if (pmb->pbval->block_bcs[outer_x1] != BoundaryFlag::block)
    ng2 = NGHOST;

  for (int n=NMASS; n<NHYDRO; ++n) {
    // left boundary
#pragma omp simd
    scale = 1.E-16;
    for (int m = -2; m <= 2; ++m) scale += w(n,k,j,il+m);
    for (int i=il; i<il+ng1; ++i) {
#pragma omp simd
      for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k,j,i+m)/(5.*scale);
      wl(n,i+1) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
      wr(n,i  ) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
      scale += w(n,k,j,i+3) - w(n,k,j,i-2);
    }

    // interior
    for (int i=il+ng1; i<=iu-ng2; ++i) {
      wl(n,i+1) = interp_cp5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i  ) = interp_cp5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
    }

    // right boundary
    scale = 1.E-16;
#pragma omp simd
    for (int m = -2; m <= 2; ++m) scale += w(n,k,j,iu-ng2+1+m);
    for (int i=iu-ng2+1; i<=iu; ++i) {
#pragma omp simd
      for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k,j,i+m)/(5.*scale);
      wl(n,i+1) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
      wr(n,i  ) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
      scale += w(n,k,j,i+3) - w(n,k,j,i-2);
    }
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X2()
//  \brief 

void Reconstruction::Weno5X2(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  MeshBlock *pmb = pmy_block_;
  Real v[5], scale;
  for (int n=0; n<NMASS; ++n) {
    scale = 1.E-16;
#pragma omp simd
    for (int m = -2; m <= 2; ++m) scale += w(n,k,j,il+m);
    for (int i=il; i<=iu; ++i) {
#pragma omp simd
      for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k,j+m,i)/(5.*scale);
      // wl(j+1/2)
      wl(n,i) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
      // wr(j-1/2)
      wr(n,i) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
      scale += w(n,k,j,i+3) - w(n,k,j,i-2);
    }
  }

  /*int ng1 = 0, ng2 = 0;
  if (pmb->pbval->block_bcs[inner_x2] != BoundaryFlag::block)
    ng1 = NGHOST;
  if (pmb->pbval->block_bcs[outer_x2] != BoundaryFlag::block)
    ng2 = NGHOST;

  if ((j < pmb->js-1 + ng1) || (j > pmb->je+1 - ng2)) {
    for (int n=NMASS; n<NHYDRO; ++n) {
      scale = 1.E-16;
#pragma omp simd
      for (int m = -2; m <= 2; ++m) scale += w(n,k,j,il+m);
      for (int i=il; i<=iu; ++i) {
#pragma omp simd
        for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k,j+m,i)/(5.*scale);
        // wl(j+1/2)
        wl(n,i) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
        // wr(j-1/2)
        wr(n,i) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
        scale += w(n,k,j,i+3) - w(n,k,j,i-2);
      }
    }
  } else*/
    for (int n=NMASS; n<NHYDRO; ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,i) = interp_cp5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
        wr(n,i) = interp_cp5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
      }
    }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X3()
//  \brief 

void Reconstruction::Weno5X3(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  MeshBlock *pmb = pmy_block_;
  Real v[5], scale;
  for (int n=0; n<NMASS; ++n) {
    scale = 1.E-16;
#pragma omp simd
    for (int m = -2; m <= 2; ++m) scale += w(n,k,j,il+m);
    for (int i=il; i<=iu; ++i) {
#pragma omp simd
      for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k+m,j,i)/(5.*scale);
      // wl(k+1/2)
      wl(n,i) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
      // wr(k-1/2)
      wr(n,i) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
      scale += w(n,k,j,i+3) - w(n,k,j,i-2);
    }
  }

  /*int ng1 = 0, ng2 = 0;
  if (pmb->pbval->block_bcs[inner_x2] != BoundaryFlag::block)
    ng1 = NGHOST;
  if (pmb->pbval->block_bcs[outer_x2] != BoundaryFlag::block)
    ng2 = NGHOST;

  if ((k < pmb->ks-1 + ng1) || (k > pmb->ke+1 - ng2)) {
    for (int n=NMASS; n<NHYDRO; ++n) {
      scale = 1.E-16;
#pragma omp simd
      for (int m = -2; m <= 2; ++m) scale += w(n,k,j,il+m);
      for (int i=il; i<=iu; ++i) {
#pragma omp simd
        for (int m = -2; m <= 2; ++m) v[m+2] = w(n,k+m,j,i)/(5.*scale);
        // wl(k+1/2)
        wl(n,i) = interp_weno5(v[4], v[3], v[2], v[1], v[0])*5.*scale;
        // wr(k-1/2)
        wr(n,i) = interp_weno5(v[0], v[1], v[2], v[3], v[4])*5.*scale;
        scale += w(n,k,j,i+3) - w(n,k,j,i-2);
      }
    }
  } else*/
    for (int n=NMASS; n<NHYDRO; ++n) {
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        wl(n,i) = interp_cp5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
        wr(n,i) = interp_cp5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
      }
    }

  return;
}
