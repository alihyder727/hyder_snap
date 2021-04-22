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

inline void ThomasTriDiag(std::vector<Eigen::Matrix<Real,5,5>>& diag, std::vector<Eigen::Matrix<Real,5,5>>& diagL,
                          std::vector<Eigen::Matrix<Real,5,5>>& diagU,std::vector<Eigen::Matrix<Real,5,1>>& rhs,
                          std::vector<Eigen::Matrix<Real,5,1>>& delta, int is, int ie)
{
   // Thomas algorithm: solve tridiagonal system
   // first row, i=is
   diag[is] = diag[is].inverse().eval();
   delta[is] = diag[is]*rhs[is];
   diag[is] *= diagU[is];
   
   // forward
   for (int i = is+1; i <= ie; ++i) {
     diag[i] = (diag[i] - diagL[i]*diag[i-1]).inverse().eval();
     delta[i] = diag[i]*(rhs[i] - diagL[i]*delta[i-1]);
     diag[i] *= diagU[i];
   }

   // update solutions, i=ie
   for (int i = ie-1; i >= is; --i) {
     delta[i] -= diag[i]*delta[i+1];
   }
}

template<typename Derived>
inline void PeriodicLinearSys(std::vector<Eigen::Matrix<Real,5,5>>& diag, std::vector<Eigen::Matrix<Real,5,5>>& diagL,
                              std::vector<Eigen::Matrix<Real,5,5>>& diagU,std::vector<Eigen::Matrix<Real,5,1>>& rhs,
                              std::vector<Eigen::Matrix<Real,5,1>>& delta,Eigen::DenseBase<Derived>& LCorner,
                              Eigen::DenseBase<Derived>& UCorner, int is, int ie)
{  // Solver for periodic linear system, ref: El-Mikkawy (2005)
   int ncells1 = ie - is + 1 + 2*NGHOST;
   // define some vectors for the solver.
   std::vector<Eigen::Matrix<Real,5,5>> alpha(ncells1);
   std::vector<Eigen::Matrix<Real,5,5>> beta(ncells1);
   std::vector<Eigen::Matrix<Real,5,5>> gamma(ncells1);
   std::vector<Eigen::Matrix<Real,5,1>> zeta(ncells1);
   std::vector<Eigen::Matrix<Real,5,5>> alpha_inv(ncells1);

   // start periodic linear system solver, ref: El-Mikkawy (2005).
   // step 1: compute alpha, beta, gamma, delta from is to ie-1.
   //                                                 -1
   // alpha[1] = a[1], gamma[1] = U, beta[1] = L*alpha[1]
   alpha[is] = diag[is]; gamma[is] = UCorner; beta[is] = LCorner;
   beta[is] *= alpha[is].inverse().eval();
   alpha_inv[is] = alpha[is].inverse().eval();

   for (int i = is+1; i < ie; i++) {
     //                             -1
     // alpha[i] = a[i] - b[i]*alpha[i-1]*c[i-1]
     alpha[i] = diag[i] - diagL[i]*alpha_inv[i-1]*diagU[i-1];
     alpha_inv[i] = alpha[i].inverse().eval();
     if (i < ie-1) {
       gamma[i] = -diagL[i]*alpha_inv[i-1]*gamma[i-1];
       beta[i] = -diagU[i-1]*alpha_inv[i]*beta[i-1];
     }
   }

   alpha[ie] = diag[ie];
   alpha_inv[ie] = alpha[ie].inverse().eval();

   for (int i = is; i < ie; i++) {
     alpha[ie] -= beta[i]*gamma[i];
   }
   gamma[ie-1] = diagU[ie-1] - diagL[ie-1]*alpha_inv[ie-2]*gamma[ie-2];
   beta[ie-1] = ( diagL[ie] - beta[ie-2]*diagU[ie-2] )*alpha_inv[ie-1];

   zeta[is] = rhs[is];
   for (int i = is+1; i < ie; i++) {
     zeta[i] = rhs[i] - diagL[i]*alpha_inv[i-1]*zeta[i-1];
   }
   zeta[ie] = rhs[ie];
   for (int i = is; i < ie; i++) {
     zeta[ie] -= beta[i]*rhs[i];
   }

   delta[ie]   = alpha_inv[ie]*zeta[ie];
   delta[ie-1] = alpha_inv[ie-1]*( zeta[ie-1] - gamma[ie-1]*delta[ie] );
   for (int i = ie-2; i >= is; --i) {
     delta[i] = alpha_inv[ie]*( zeta[i] - diagU[i]*delta[i+1] - gamma[i]*delta[ie] );
   }

}

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
        PeriodicForwardSweep(a, b, c, delta, corr, dt, k, j, is, ie);
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

  delete [] gamma_m1;
}
