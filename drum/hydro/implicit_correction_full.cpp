//! \file implicit_correction.cpp
//  \brief vertical implicit roe solver

// C/C++ headers
#include <iostream>
#include <vector>

// Eigen headers
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/Dense"

// Athena++ headers
#include "../math/core.h"  // _sqr
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "hydro.hpp"

inline void RoeAverage(Real prim[], Real gm1, AthenaArray<Real> const& w, int k, int j, int i)
{
  Real wl[NHYDRO], wr[NHYDRO];
  for (int n = 0; n < NHYDRO; ++n) {
    wl[n] = w(n,k,j,i-1);
    wr[n] = w(n,k,j,i);
  }
  Real sqrtdl = sqrt(wl[IDN]);
  Real sqrtdr = sqrt(wr[IDN]);
  Real isdlpdr = 1.0/(sqrtdl + sqrtdr);

  // Roe average scheme
  // Flux in the interface between i-th and i+1-th cells: 
  // A(i+1/2) = [sqrt(rho(i))*A(i) + sqrt(rho(i+1))*A(i+1)]/(sqrt(rho(i)) + sqrt(rho(i+1)))

  prim[IDN] = sqrtdl*sqrtdr;					
  prim[IVX] = (sqrtdl*wl[IVX] + sqrtdr*wr[IVX])*isdlpdr;
  prim[IVY] = (sqrtdl*wl[IVY] + sqrtdr*wr[IVY])*isdlpdr;
  prim[IVZ] = (sqrtdl*wl[IVZ] + sqrtdr*wr[IVZ])*isdlpdr;

  for (int i = 1; i < NMASS; ++i)
    prim[i] = (sqrtdl*wl[i] + sqrtdr*wr[i])*isdlpdr;

  // Etot of the left side.
  Real el = wl[IPR]/gm1 + 0.5*wl[IDN]*(sqr(wl[IVX]) + sqr(wl[IVY]) + sqr(wl[IVZ]));

  // Etot of the right side.
  Real er = wr[IPR]/gm1 + 0.5*wr[IDN]*(sqr(wr[IVX]) + sqr(wr[IVY]) + sqr(wr[IVZ]));

  // Enthalpy divided by the density.
  Real hbar = ((el + wl[IPR])/sqrtdl + (er + wr[IPR])/sqrtdr)*isdlpdr;

  // Roe averaged pressure
  prim[IPR] = (hbar - 0.5*(sqr(prim[IVX]) + sqr(prim[IVY]) + sqr(prim[IVZ])))
              *gm1/(gm1 + 1.)*prim[IDN];
}

template<typename Derived>
inline void Eigenvalue(Eigen::MatrixBase<Derived>& Lambda, Real u, Real cs)
{
  Lambda << u-cs,  0.,    0.,    0.,    0.,
            0.,    u,     0.,    0.,    0.,
            0.,    0.,    u+cs,  0.,    0.,
            0.,    0.,    0.,    u,     0.,
            0.,    0.,    0.,    0.,    u;

  Lambda = Lambda.cwiseAbs();
}

template<typename Derived1, typename Derived2>
inline void Eigenvector(Eigen::DenseBase<Derived1>& Rmat, Eigen::DenseBase<Derived2>& Rimat,
  Real prim[], Real cs, Real gm1)
{
  Real r = prim[IDN];
  Real u = prim[IVX];
  Real v = prim[IVY];
  Real w = prim[IVZ];
  Real p = prim[IPR];
  
  Real ke = 0.5*(u*u + v*v + w*w);
  Real hp = (gm1 + 1.)/gm1*p/r;
  Real h = hp + ke;
  
  Rmat << 1.,     1.,     1.,      0.,     0.,
          u-cs,   u,      u+cs,    0.,     0.,
          v,      v,      v,       1.,     0.,
          w,      w,      w,       0.,     1.,
          h-u*cs, ke,     h+u*cs,  v,      w;

  Rimat << (cs*ke+u*hp)/(2.*cs*hp),  (-hp-cs*u)/(2.*cs*hp),  -v/(2.*hp),  -w/(2.*hp),  1./(2.*hp),
           (hp-ke)/hp,               u/hp,                   v/hp,        w/hp,        -1./hp,
           (cs*ke-u*hp)/(2.*cs*hp),  (hp-cs*u)/(2.*cs*hp),   -v/(2.*hp),  -w/(2.*hp),  1./(2.*hp),
           -v,                       0.,                     1.,          0.,          0.,
           -w,                       0.,                     0.,          1.,          0.;
}

template<typename Derived1>
inline void FluxJacobian(Eigen::DenseBase<Derived1>& dfdq,
  Real gm1, AthenaArray<Real> const& w, int k, int j, int i)
{
  // flux derivative
  // Input variables are density, velocity field and energy.
  // The primitives of cell (n,i)
  Real v1  = w(IVX,k,j,i);
  Real v2  = w(IVY,k,j,i);
  Real v3  = w(IVZ,k,j,i); 
  Real rho = w(IDN,k,j,i);
  Real pres = w(IPR,k,j,i);
  Real s2 = v1*v1 + v2*v2 + v3*v3;

  Real c1 = ((gm1-1)*s2/2 - (gm1+1)/gm1*pres/rho)*v1;
  Real c2 = (gm1+1)/gm1*pres/rho + s2/2 - gm1*v1*v1;

  dfdq  << 0,                 1.,               0.,           0.,           0.,
           gm1*s2/2-v1*v1,    (2.-gm1)*v1,      -gm1*v2,      -gm1*v3,      gm1,
           -v1*v2,            v2,               v1,           0.,           0.,
           -v1*v3,            v3,               0.,           v1,           0.,
           c1,                c2,               -gm1*v2*v1,   -gm1*v3*v1,   (gm1+1)*v1;
}

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

void Hydro::ImplicitCorrectionFull(AthenaArray<Real> &du, AthenaArray<Real> const& w, Real dt)
{ 
  // Implicit Correction for i-direction.
  MeshBlock *pmb = pmy_block;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  int ks = pmb->ks; int js = pmb->js; int is = pmb->is;
  int ke = pmb->ke; int je = pmb->je; int ie = pmb->ie;

  int ncells1 = ie - is + 1 + 2*NGHOST;

  AthenaArray<Real> &vol = cell_volume_;
  AthenaArray<Real> &farea = x1face_area_;

  // eigenvectors, eigenvalues, inverse matrix of eigenvectors.
  Eigen::Matrix<Real,5,5> Rmat, Lambda, Rimat;

  // reduced diffusion matrix |A_{i-1/2}|, |A_{i+1/2}|
  Eigen::Matrix<Real,5,5> Am, Ap;

  Real prim[NHYDRO]; // Roe averaged primitive variables of cell i-1/2

  std::vector<Eigen::Matrix<Real,5,1>> rhs(ncells1);
  std::vector<Eigen::Matrix<Real,5,5>> a(ncells1), b(ncells1), c(ncells1);
  std::vector<Eigen::Matrix<Real,5,1>> delta(ncells1);
  std::vector<Eigen::Matrix<Real,5,5>> dfdq(ncells1);

  // 0. forcing and volume matrix
  Real gamma = pmb->peos->GetGamma();
  Real grav = hsrc.GetG1();
  Eigen::Matrix<Real,5,5> Phi, Dt, Bnd;
  Phi  << 0.,    0.,    0.,    0.,    0.,
          grav,  0.,    0.,    0.,    0., 
          0.,    0.,    0.,    0.,    0.,
          0.,    0.,    0.,    0.,    0., 
          0.,    grav,  0.,    0.,    0.;

  Dt   << 1./dt, 0.,    0.,    0.,    0.,
          0., 1./dt,    0.,    0.,    0.,
          0.,    0., 1./dt,    0.,    0.,
          0.,    0.,    0., 1./dt,    0.,
          0.,    0.,    0.,    0., 1./dt;

  Bnd  << 1.,    0.,    0.,    0.,    0.,
          0.,   -1.,    0.,    0.,    0.,
          0.,    0.,    1.,    0.,    0.,
          0.,    0.,    0.,    1.,    0.,
          0.,    0.,    0.,    0.,    1.;

  Real *gamma_m1 = new Real [ncells1];

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->Face1Area(k, j, is, ie+1, farea);
      pcoord->CellVolume(k, j, is, ie, vol);

      // 3. calculate and save flux Jacobian matrix
      for (int i = is-1; i <= ie+1; ++i) {
        //FluxJacobian(dfdq[i], gamma - 1, w, k, j, i);
        Real fsig = 1., feps = 1.;
        for (int n = 1 + NVAPOR; n < NMASS; ++n) {
          fsig += w(n,k,j,i)*(pthermo->GetCvRatio(n) - 1.);
          feps -= w(n,k,j,i);
        }
        for (int n = 1; n <= NVAPOR; ++n) {
          fsig += w(n,k,j,i)*(pthermo->GetCvRatio(n) - 1.);
          feps += w(n,k,j,i)*(1./pthermo->GetMassRatio(n) - 1.);
        }

        gamma_m1[i] = (gamma - 1.)*feps/fsig;
        FluxJacobian(dfdq[i], gamma_m1[i], w, k, j, i);
      }

      // 5. set up diffusion matrix and tridiagonal coefficients
      // left edge
      Real gm1 = 0.5*(gamma_m1[is-1] + gamma_m1[is]);
      RoeAverage(prim, gm1, w, k, j, is);
      Real cs = pmb->peos->SoundSpeed(prim);
      Eigenvalue(Lambda, prim[IVX], cs);
      Eigenvector(Rmat, Rimat, prim, cs, gm1);
      Am = Rmat*Lambda*Rimat;

      for (int i = is; i <= ie; ++i) {
        // right edge
        gm1 = 0.5*(gamma_m1[i] + gamma_m1[i+1]);
        RoeAverage(prim, gm1, w, k, j, i+1);
        Real cs = pmb->peos->SoundSpeed(prim);
        Eigenvector(Rmat, Rimat, prim, cs, gm1);
        Ap = Rmat*Lambda*Rimat;

        // set up diagonals a, b, c.
        a[i] = (Am*farea(i) + Ap*farea(i+1) + (farea(i+1) - farea(i))*dfdq[i])/(2.*vol(i)) 
               + Dt - Phi;
        b[i] = -(Am + dfdq[i-1])*farea(i)/(2.*vol(i));
        c[i] = -(Ap - dfdq[i+1])*farea(i+1)/(2.*vol(i));

        // Shift one cell: i -> i+1
        Am = Ap;
      }

      // 5. fix boundary condition
      if (pmb->pbval->block_bcs[inner_x1] != BoundaryFlag::periodic) {
        a[is] += b[is]*Bnd;
      } 
      if (pmb->pbval->block_bcs[outer_x1] != BoundaryFlag::periodic) {
        a[ie] += c[ie]*Bnd;
      }

      // 6. solve tridiagonal system
      for (int i = is; i <= ie; ++i) {
        rhs[i](0) = du(IDN,k,j,i)/dt;
        for (int n = NMASS; n <= NMASS+3; ++n) {
          rhs[i](n-NMASS+1) = du(n,k,j,i)/dt;
        }
      }

      if ( (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::periodic) &&
           (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::periodic) ) {
        // LU decomposition for periodic (i.e., circuit) linear system
        PeriodicLinearSys(a,b,c,rhs,delta,c[ie],b[is],is,ie);
      } else {
        // TDMA for tridiagonal linear system
        ThomasTriDiag(a,b,c,rhs,delta,is,ie);
      }

      // 7. update conserved variables, i = ie
      for (int i = is; i <= ie; ++i) {
        //for (int n = 0; n < 5; ++n)
        //  du(n,k,j,i) = delta[i](n);
        du(IDN,k,j,i) = delta[i](0);
        for (int n = NMASS; n <= NMASS+3; ++n) {
          du(n,k,j,i) = delta[i](n-NMASS+1);
        }
      }
    }

  delete [] gamma_m1;
}
