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

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::Weno5X1()
//  \brief 

void Reconstruction::Weno5X1(const int k, const int j, const int il, const int iu,
  const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
  AthenaArray<Real> &wl, AthenaArray<Real> &wr)
{
  for (int n=0; n<=NVAPOR; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i+1) = interp_weno5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
      wr(n,i) = interp_weno5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
    }

  for (int n=NVAPOR+1; n<NMASS; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i+1) = interp_weno3(w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1));
      wr(n,i) = interp_weno3(w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1));
    }

  for (int n=NMASS; n<NHYDRO; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i+1) = interp_cp5(w(n,k,j,i+2),w(n,k,j,i+1),w(n,k,j,i),w(n,k,j,i-1),w(n,k,j,i-2));
      wr(n,i) = interp_cp5(w(n,k,j,i-2),w(n,k,j,i-1),w(n,k,j,i),w(n,k,j,i+1),w(n,k,j,i+2));
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
  for (int n=0; n<=NVAPOR; ++n)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = interp_weno5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
      wr(n,i) = interp_weno5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
    }

  for (int n=NVAPOR+1; n<NMASS; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = interp_weno3(w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i));
      wr(n,i) = interp_weno3(w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i));
    }

  for (int n=NMASS; n<NHYDRO; ++n)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = interp_cp5(w(n,k,j+2,i),w(n,k,j+1,i),w(n,k,j,i),w(n,k,j-1,i),w(n,k,j-2,i));
      wr(n,i) = interp_cp5(w(n,k,j-2,i),w(n,k,j-1,i),w(n,k,j,i),w(n,k,j+1,i),w(n,k,j+2,i));
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
  for (int n=0; n<=NVAPOR; ++n)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = interp_weno5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
      wr(n,i) = interp_weno5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
    }

  for (int n=NVAPOR+1; n<NMASS; ++n) 
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = interp_weno3(w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i));
      wr(n,i) = interp_weno3(w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i));
    }

  for (int n=NMASS; n<NHYDRO; ++n)
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      wl(n,i) = interp_cp5(w(n,k+2,j,i),w(n,k+1,j,i),w(n,k,j,i),w(n,k-1,j,i),w(n,k-2,j,i));
      wr(n,i) = interp_cp5(w(n,k-2,j,i),w(n,k-1,j,i),w(n,k,j,i),w(n,k+1,j,i),w(n,k+2,j,i));
    }

  return;
}
