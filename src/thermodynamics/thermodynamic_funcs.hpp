#ifndef THERMODYNAMIC_FUNCS_HPP
#define THERMODYNAMIC_FUNCS_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "../athena.hpp"

// frequently used inline functions
inline Real SatVaporPresIdeal(Real t, Real p3, Real beta, Real delta) {
  return p3*exp((1. - 1./t)*beta - delta*log(t));
}

// alpha = L/cv evaluated at current temperature
Real GasCloudIdeal(Real const q[], int iv, int ic,
  Real t3, Real p3, Real alpha, Real beta, Real delta);

Real q_gas(Real const w[]);
Real q_eps(Real const w[], Real const eps[]);
Real q_sig(Real const w[], Real const sig[]);
Real qhat_eps(Real const q[], Real const eps[]);
//Real qhat_rcp(Real const q[], Real const rcp[]);

void mass2molar(Real q[], Real const w[], Real const eps[], Real Rd);
void molar2mass(Real w[], Real const q[], Real const eps[], Real Rd);

Real dlnTdlnP(Real const q[], int const isat[], 
  Real const rcp[], Real const beta[], Real const delta[], Real const t3[], Real gamma);

void rk4_integrate_z(Real q[], int isat[], Real rcp[], Real const eps[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma,
  Real g_ov_Rd, Real dz, int method, Real user_dTdz = 0);

void rk4_integrate_z_adaptive(
  Real q[], int isat[], Real rcp[], Real const eps[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma,
  Real g_ov_Rd, Real dz, Real ftol, int method, Real user_dTdz = 0);

void rk4_integrate_lnp(Real prim[], int isat[], Real const rcp[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma, Real dp);

void rk4_integrate_lnp_adaptive(Real prim[], int isat[], Real const rcp[],
  Real const beta[], Real const delta[], Real const t3[], Real const p3[], Real gamma, Real dp, Real ftol);

void update_gamma(Real& gamma, Real rcp[], Real const prim[]);

bool read_idealized_parameters(std::string fname, Real &pmin, Real &pmax, AthenaArray<Real> &dTdz, 
  Real* &pmin_q, Real* &pmax_q, AthenaArray<Real>* &dqdz, std::vector<int> &mindex);

#endif
