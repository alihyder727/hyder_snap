#ifndef FLUX_DECOMPOSITION_HPP
#define FLUX_DECOMPOSITION_HPP

// Athena++ headers
#include "../../math/core.h"  // sqr
#include "../../athena.hpp"

// Eigen headers
#include "../../math/eigen335/Eigen/Core"
#include "../../math/eigen335/Eigen/Dense"

inline void RoeAverage(Real prim[], Real gm1, Real wl[], Real wr[])
{
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
inline void FluxJacobian(Eigen::DenseBase<Derived1>& dfdq, Real gm1, Real w[])
{
  // flux derivative
  // Input variables are density, velocity field and energy.
  // The primitives of cell (n,i)
  Real v1  = w[IVX];
  Real v2  = w[IVY];
  Real v3  = w[IVZ];
  Real rho = w[IDN];
  Real pres = w[IPR];
  Real s2 = v1*v1 + v2*v2 + v3*v3;

  Real c1 = ((gm1-1)*s2/2 - (gm1+1)/gm1*pres/rho)*v1;
  Real c2 = (gm1+1)/gm1*pres/rho + s2/2 - gm1*v1*v1;

  dfdq  << 0,                 1.,               0.,           0.,           0.,
           gm1*s2/2-v1*v1,    (2.-gm1)*v1,      -gm1*v2,      -gm1*v3,      gm1,
           -v1*v2,            v2,               v1,           0.,           0.,
           -v1*v3,            v3,               0.,           v1,           0.,
           c1,                c2,               -gm1*v2*v1,   -gm1*v3*v1,   (gm1+1)*v1;
}

#endif
