//! \file implicit_correction.cpp
//  \brief vertical implicit roe solver

// C/C++ headers
#include <iostream>
#include <vector>

// Eigen headers
#include "../math/eigen335/Eigen/Core"
#include "../math/eigen335/Eigen/Dense"

// Athena++ headers
#include "../math/core.h"  // sqr
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "vertical_communication.hpp"

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

template<typename Derived1, typename Derived2>
inline void FluxJacobian(Eigen::DenseBase<Derived1>& dfdq1,
                         Eigen::DenseBase<Derived2>& dfdq2,
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

  dfdq1 <<  0.,         0.,
            -gm1*v2,    -gm1*v3,
            -gm1*v2*v1, -gm1*v3*v1;

  dfdq2 <<  0,                1.,           0., 
            gm1*s2/2-v1*v1/2, (2.-gm1)*v1,  gm1,
            c1,               c2,           (gm1+1)*v1;

}

void Hydro::ImplicitCorrectionReduced(AthenaArray<Real> &du, AthenaArray<Real> const& w, Real dt)
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
  Eigen::Matrix<Real,3,2> Am1, Ap1;
  Eigen::Matrix<Real,3,3> Am2, Ap2;
  Eigen::Matrix<Real,3,1> rhs, sm, sp;
  Eigen::Matrix<Real,2,1> dqm, dqp;

  Real prim[NHYDRO]; // Roe averaged primitive variables of cell i-1/2

  std::vector<Eigen::Matrix<Real,3,3>> a(ncells1), b(ncells1), c(ncells1);
  std::vector<Eigen::Matrix<Real,3,1>> delta(ncells1);

  std::vector<Eigen::Matrix<Real,3,2>> dfdq1(ncells1);
  std::vector<Eigen::Matrix<Real,3,3>> dfdq2(ncells1);

  std::vector<Eigen::Matrix<Real,3,1>> corr(ncells1);

  // 0. forcing and volume matrix
  Real gamma = pmb->peos->GetGamma();
  Real grav = hsrc.GetG1();
  Eigen::Matrix<Real,3,3> Phi, Dt, Bnd;
  Phi  << 0.,    0.,    0.,
          grav,  0.,    0.,
          0.,    grav,  0.;

  Dt   << 1./dt, 0.,    0.,
          0., 1./dt,    0.,
          0.,    0., 1./dt;

  Bnd  << 1.,    0.,    0.,
          0.,   -1.,    0.,
          0.,    0.,    1.;

  Real *gamma_m1 = new Real [ncells1];
  int idn = 0, ivx = 1, ivy = 2, ivz = 3, ien = 4;

  pvc->FindNeighbors();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->Face1Area(k, j, is, ie+1, farea);
      pcoord->CellVolume(k, j, is, ie, vol);

      // 3. calculate and save flux Jacobian matrix
      for (int i = is-1; i <= ie+1; ++i) {
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
        FluxJacobian(dfdq1[i], dfdq2[i], gamma_m1[i], w, k, j, i);
      }

      // 5. set up diffusion matrix and tridiagonal coefficients
      // left edge
      Real gm1 = 0.5*(gamma_m1[is-1] + gamma_m1[is]);
      RoeAverage(prim, gm1, w, k, j, is);
      Real cs = pmb->peos->SoundSpeed(prim);
      Eigenvalue(Lambda, prim[IVX], cs);
      Eigenvector(Rmat, Rimat, prim, cs, gm1);
      Am = Rmat*Lambda*Rimat;
      Am1 << Am(idn,ivy), Am(idn,ivz),
             Am(ivx,ivy), Am(ivx,ivz),
             Am(ien,ivy), Am(ien,ivz);
      Am2 << Am(idn,idn), Am(idn,ivx), Am(idn,ien),
             Am(ivx,idn), Am(ivx,ivx), Am(ivx,ien),
             Am(ien,idn), Am(ien,ivx), Am(ien,ien);

      for (int i = is; i <= ie; ++i) {
        // right edge
        gm1 = 0.5*(gamma_m1[i] + gamma_m1[i+1]);
        RoeAverage(prim, gm1, w, k, j, i+1);
        Real cs = pmb->peos->SoundSpeed(prim);
        Eigenvector(Rmat, Rimat, prim, cs, gm1);
        Ap = Rmat*Lambda*Rimat;
        Ap1 << Ap(idn,ivy), Ap(idn,ivz),
               Ap(ivx,ivy), Ap(ivx,ivz),
               Ap(ien,ivy), Ap(ien,ivz);
        Ap2 << Ap(idn,idn), Ap(idn,ivx), Ap(idn,ien),
               Ap(ivx,idn), Ap(ivx,ivx), Ap(ivx,ien),
               Ap(ien,idn), Ap(ien,ivx), Ap(ien,ien);

        // set up diagonals a, b, c.
        a[i] = (Am2*farea(i) + Ap2*farea(i+1) + (farea(i+1) - farea(i))*dfdq2[i])/(2.*vol(i)) 
               + Dt - Phi;
        b[i] = -(Am2 + dfdq2[i-1])*farea(i)/(2.*vol(i));
        c[i] = -(Ap2 - dfdq2[i+1])*farea(i+1)/(2.*vol(i));

        // flux correction
        dqm << du(IVY,k,j,i), du(IVZ,k,j,i);
        dqp << du(IVY,k,j,i+1), du(IVZ,k,j,i+1);
        sm = 0.5*((dfdq1[i-1] + Am1)*dqm + (dfdq1[i] - Am1)*dqp);
        sp = 0.5*((dfdq1[i] + Ap1)*dqm + (dfdq1[i+1] - Ap1)*dqp);
        corr[i] = (sp*farea(i+1) - sm*farea(i))/vol(i);

        // Shift one cell: i -> i+1
        Am1 = Ap1;
        Am2 = Ap2;
      }

      // 5. fix boundary condition
      if (!pvc->has_bot_neighbor)
        a[is] += b[is]*Bnd;
      if (!pvc->has_top_neighbor)
        a[ie] += c[ie]*Bnd;

      // 6. Thomas algorithm: solve tridiagonal system
      // first row, i=is
      rhs(0) = du(IDN,k,j,is)/dt;
      rhs(1) = du(IVX,k,j,is)/dt;
      rhs(2) = du(IEN,k,j,is)/dt;
      rhs -= corr[is];

      if (!pvc->has_bot_neighbor) {
        a[is] = a[is].inverse().eval();
        delta[is] = a[is]*rhs;
        a[is] *= c[is];
      } else {
        pvc->RecvBotImplicitData(a[is-1], delta[is-1], k, j, pvc->bblock);
        a[is] = (a[is] - b[is]*a[is-1]).inverse().eval();
        delta[is] = a[is]*(rhs - b[is]*delta[is-1]);
        a[is] *= c[is];
      }
       
      // forward
      for (int i = is+1; i <= ie; ++i) {
        rhs(0) = du(IDN,k,j,i)/dt;
        rhs(1) = du(IVX,k,j,i)/dt;
        rhs(2) = du(IEN,k,j,i)/dt;
        rhs -= corr[i];

        a[i] = (a[i] - b[i]*a[i-1]).inverse().eval();
        delta[i] = a[i]*(rhs - b[i]*delta[i-1]);
        a[i] *= c[i];
      }

      pvc->SaveImplicitData(a, delta, k, j, is, ie);

      if (pvc->has_top_neighbor)
        pvc->SendTopImplicitData(a[ie], delta[ie], k, j, pvc->tblock);
    }

  // back substitution
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pvc->LoadImplicitData(a, delta, k, j, is, ie);
      if (pvc->has_top_neighbor) {
        pvc->RecvTopImplicitData(delta[ie+1], k, j, pvc->tblock);
        delta[ie] -= a[ie]*delta[ie+1];
      }

      // update solutions, i=ie
      for (int i = ie-1; i >= is; --i)
        delta[i] -= a[i]*delta[i+1];

      // 7. update conserved variables, i = ie
      for (int i = is; i <= ie; ++i) {
        du(IDN,k,j,i) = delta[i](0);
        du(IVX,k,j,i) = delta[i](1);
        du(IEN,k,j,i) = delta[i](2);
      }

      if (pvc->has_bot_neighbor)
        pvc->SendBotImplicitData(delta[is], k, j, pvc->bblock);
    }

  if (pvc->has_top_neighbor)
    pvc->WaitSendTopImplicit(ks, ke, js, je);

  if (pvc->has_bot_neighbor)
    pvc->WaitSendBotImplicit(ks, ke, js, je);

  delete [] gamma_m1;
}
