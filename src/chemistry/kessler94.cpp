// Athena++ header
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "chemistry.hpp"

// Implementation of Kessler 1994 microphysics scheme 
// For each condensable species, the first phase is gas; the second phase is cloud
// and the third phase is precipitation
// Three reactions:
// 1. autoconversion
// 2. accretion
// 3. evaporation

// Because the mixing ratios may turn negative for a variety of reasons
// Make sure that the chemical scheme works with negative mixing ratios
void Kessler94::AssembleReactionMatrix(Real *r0, Real **r1, Real const q[], Real time)
{
  /*Real k1 = kc_.at(1);
  Real k2 = kc_.at(2);
  Real k3 = kc_.at(3);
  std::fill(*r1_, *r1_ + NMASS*NMASS, 0.);
  std::fill(r0_, r0_ + NMASS, 0.);

  Real qtol = 1.; // total gas mols
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    qtol -= q[n];

  for (int n = 1; n <= NVAPOR; ++n) {
    int iv = n, ic = n + NVAPOR, ip = n + 2*NVAPOR;
    Real svp = pmy_block_->pthermo->SatVaporPressure(q[IDN], ic);
    Real dq = q[n] - svp/q[IPR]*qtol; // q - qs

    if (dq < 0.) {  // evaporation
      r0[iv] += -k3*q[ip]*dq;
      r0[ip] += k3*q[ip]*dq;
      r1[iv][iv] += -k3*q[ip];
      r1[iv][ip] += -k3*dq;
      r1[ip][iv] += k3*q[ip];
      r1[ip][ip] += k3*dq;
    }

    // autoconversion
    r0[ic] += -k1*q[ic];
    r0[ip] += k1*q[ic];
    r1[ic][ic] += -k1;
    r1[ip][ic] += k1;

    // accretion
    r0[ic] += -k2*q[ic]*q[ip];
    r0[ip] += k2*q[ic]*q[ip];
    r1[ic][ic] += -k2*q[ip];
    r1[ic][ip] += -k2*q[ic];
    r1[ip][ic] += k2*q[ip];
    r1[ip][ip] += k2*q[ic];
  }*/
}
