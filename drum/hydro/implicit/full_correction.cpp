//! \file full_correction.cpp
//  \brief vertical implicit roe solver

// C/C++ headers
#include <iostream>
#include <vector>

// Eigen headers
#include "../../math/eigen335/Eigen/Core"
#include "../../math/eigen335/Eigen/Dense"

// Athena++ headers
#include "../../math/core.h"  // _sqr
#include "../../mesh/mesh.hpp"
#include "../../eos/eos.hpp"
#include "../../thermodynamics/thermodynamics.hpp"
#include "../hydro.hpp"
#include "flux_decomposition.hpp"
#include "implicit_solver.hpp"
#include "forward_backward.hpp"
#include "periodic_forward_backward.hpp"
#include "../jacobian_functions.hpp"


void ImplicitSolver::FullCorrection(AthenaArray<Real>& du,
  AthenaArray<Real> const& w, Real dt)
{ 
  MeshBlock *pmb = pmy_hydro->pmy_block;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  int is, ie, js, je, ks, ke;
  int idn = 0, ivx = 1, ivy = 2, ivz = 3, ien = 4;
  if (mydir == X1DIR) {
    ks = pmb->ks, js = pmb->js, is = pmb->is;
    ke = pmb->ke, je = pmb->je, ie = pmb->ie;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i)
            du_(n,k,j,i) = du(n,k,j,i);
  } else if (mydir == X2DIR) {
    ks = pmb->is, js = pmb->ks, is = pmb->js;
    ke = pmb->ie, je = pmb->ke, ie = pmb->je;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i)
            du_(n,k,j,i) = du(n,j,i,k);
  } else { // X3DIR
    ks = pmb->js, js = pmb->is, is = pmb->ks;
    ke = pmb->je, je = pmb->ie, ie = pmb->ke;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          for (int i = is; i <= ie; ++i)
            du_(n,k,j,i) = du(n,i,k,j);
  }

  // eigenvectors, eigenvalues, inverse matrix of eigenvectors.
  Eigen::Matrix<Real,5,5> Rmat, Lambda, Rimat;

  // reduced diffusion matrix |A_{i-1/2}|, |A_{i+1/2}|
  Eigen::Matrix<Real,5,5> Am, Ap;

  Real prim[NHYDRO]; // Roe averaged primitive variables of cell i-1/2

  int nc = ie - is + 1 + 2*NGHOST;
  std::vector<Eigen::Matrix<Real,5,1>> rhs(nc);
  std::vector<Eigen::Matrix<Real,5,5>> a(nc), b(nc), c(nc);
  std::vector<Eigen::Matrix<Real,5,1>> delta(nc);
  std::vector<Eigen::Matrix<Real,5,5>> dfdq(nc);
  std::vector<Eigen::Matrix<Real,5,1>> corr(nc);  // place holder

  // 0. forcing and volume matrix
  FindNeighbors();

  Real gamma = pmb->peos->GetGamma();
  Eigen::Matrix<Real,5,5> Phi, Dt, Bnds, Bnde;

  Dt.setIdentity();
  Dt *= 1./dt;

  Bnds.setIdentity();
  Bnds(ivx,ivx) = -1;
  if (pole_at_bot)
    Bnds(ivy,ivy) = -1;

  Bnde.setIdentity();
  Bnde(ivx,ivx) = -1;
  if (pole_at_top)
    Bnde(ivy,ivy) = -1;

  Real *gamma_m1 = new Real [nc];

  Real wl[NHYDRO], wr[NHYDRO];
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      // 3. calculate and save flux Jacobian matrix
      for (int i = is-1; i <= ie+1; ++i) {
        Real fsig = 1., feps = 1.;
        CopyPrimitives(wl, wr, w, k, j, i, mydir);
        for (int n = 1 + NVAPOR; n < NMASS; ++n) {
          fsig += w(n,k,j,i)*(pthermo->GetCvRatio(n) - 1.);
          feps -= w(n,k,j,i);
        }
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += w(n,k,j,i)*(pthermo->GetCvRatio(n) - 1.);
          feps += w(n,k,j,i)*(1./pthermo->GetMassRatio(n) - 1.);
        }

        gamma_m1[i] = (gamma - 1.)*feps/fsig;
        FluxJacobian(dfdq[i], gamma_m1[i], wr, mydir);
      } // 5. set up diffusion matrix and tridiagonal coefficients
      // left edge
      CopyPrimitives(wl, wr, w, k, j, is, mydir);
      Real gm1 = 0.5*(gamma_m1[is-1] + gamma_m1[is]);
      RoeAverage(prim, gm1, wl, wr);
      Real cs = pmb->peos->SoundSpeed(prim);
      Eigenvalue(Lambda, prim[IVX+mydir], cs);
      Eigenvector(Rmat, Rimat, prim, cs, gm1, mydir);
      Am = Rmat*Lambda*Rimat;

      for (int i = is-1; i <= ie+1; ++i) {
        CopyPrimitives(wl, wr, w, k, j, i+1, mydir);
        // right edge
        gm1 = 0.5*(gamma_m1[i] + gamma_m1[i+1]);
        RoeAverage(prim, gm1, wl, wr);
        Real cs = pmb->peos->SoundSpeed(prim);
        Eigenvalue(Lambda, prim[IVX+mydir], cs);
        Eigenvector(Rmat, Rimat, prim, cs, gm1, mydir);
        Ap = Rmat*Lambda*Rimat;

        // set up diagonals a, b, c, and Jacobian of the forcing function
        Real aleft, aright, vol;
        if (mydir == X1DIR) {
          aleft = pcoord->GetFace1Area(k,j,i);
          aright = pcoord->GetFace1Area(k,j,i+1);
          vol = pcoord->GetCellVolume(k,j,i);
          JACOBIAN_FUNCTION(Phi,k,j,i,pmy_hydro,mydir);
        } else if (mydir == X2DIR) {
          aleft = pcoord->GetFace2Area(j,i,k);
          aright = pcoord->GetFace2Area(j,i+1,k);
          vol = pcoord->GetCellVolume(j,i,k);
          JACOBIAN_FUNCTION(Phi,j,i,k,pmy_hydro,mydir);
        } else { // X3DIR
          aleft = pcoord->GetFace3Area(i,k,j);
          aright = pcoord->GetFace3Area(i+1,k,j);
          vol = pcoord->GetCellVolume(i,k,j);
          JACOBIAN_FUNCTION(Phi,i,k,j,pmy_hydro,mydir);
        }
        SaveForcingJacobian(Phi,k,j,i);

        a[i] = (Am*aleft + Ap*aright + (aright - aleft)*dfdq[i])/(2.*vol) 
               + Dt - Phi;
        b[i] = -(Am + dfdq[i-1])*aleft/(2.*vol);
        c[i] = -(Ap - dfdq[i+1])*aright/(2.*vol);

        // Shift one cell: i -> i+1
        Am = Ap;
      }

      // 5. fix boundary condition
      if (first_block && !periodic_boundary)
        a[is] += b[is]*Bnds;

      if (last_block && !periodic_boundary)
        a[ie] += c[ie]*Bnde;

      // 6. solve tridiagonal system
      if (periodic_boundary)
        PeriodicForwardSweep(a, b, c, corr, dt, k, j, is, ie);
      else
        ForwardSweep(a, b, c, delta, corr, dt, k, j, is, ie);
    }

  if (periodic_boundary)
    PeriodicBackwardSubstitution(a, c, delta, ks, ke, js, je, is, ie);
  else
    BackwardSubstitution(a, delta, ks, ke, js, je, is, ie);

  if (mydir == X1DIR) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,k,j,i) = du_(IDN,k,j,i);
          for (int n = NMASS; n <= NMASS + 3; ++n)
            du(n,k,j,i) = du_(n,k,j,i);
        }
  } else if (mydir == X2DIR) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,j,i,k) = du_(IDN,k,j,i);
          for (int n = NMASS; n <= NMASS + 3; ++n)
            du(n,j,i,k) = du_(n,k,j,i);
        }
  } else {  // X3DIR
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,i,k,j) = du_(IDN,k,j,i);
          for (int n = NMASS; n <= NMASS + 3; ++n)
            du(n,i,k,j) = du_(n,k,j,i);
        }
  }
  pmy_hydro->implicit_done = this;

  delete [] gamma_m1;
}
