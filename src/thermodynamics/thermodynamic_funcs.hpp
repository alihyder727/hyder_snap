/** @file thermodynamic_funcs.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday Jun 01, 2021 15:42:38 PDT
 * @bug No known bugs.
 */

#ifndef THERMODYNAMICS_FUNCS_HPP
#define THERMODYNAMICS_FUNCS_HPP

#include "thermodynamics.hpp"

//! Potential temperature
template<typename T>
Real PotentialTemp(T w, Real p0, Thermodynamics *pthermo) {
  Real chi = pthermo->GetChi(w);
  Real temp = pthermo->GetTemp(w);
  return temp*pow(p0/w[IPR], chi);
}

//! Moist static energy
template<typename A, typename B>
Real MoistStaticEnergy(A const w, B const c, Real gz, Thermodynamics *pthermo) {
  Real temp = pthermo->GetTemp(w);
  Real IE = w[IDN]*pthermo->GetMeanCp(w)*temp;
  Real rho = w[IDN];
  for (int n = 0; n < 2*NVAPOR; ++n) {
    rho += c[n];
    IE += c[n]*pthermo->GetCp(1+NVAPOR+n)*temp;
    IE -= c[n]*pthermo->GetLatent(1+NVAPOR+n);
  }
  return IE/rho + gz;
}

//! Relative humidity
template<typename T>
Real RelativeHumidity(T w, int iv, Thermodynamics *pthermo) {
  Real dw[1+NVAPOR];
  pthermo->SaturationSurplus(dw, w);
  return w[iv]/(w[iv] - dw[iv]);
}

//! Equivalent potential temperature
template<typename A, typename B>
Real MoistEntropy(A const, B const p0, Thermodynamics *pthermo) {
/*#if (NVAPOR > 0)
  Real gamma = pthermo->pmy_block->peos->GetGamma();
  Real tem[1] = {pthermo->GetTemp(prim)};
  update_gamma(gamma, tem);
  Real cpd = Rd_*gamma/(gamma - 1.);
  Real temp = tem[0]
  Real pres = prim[IPR];

  Real qd = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    qd -= prim[n];

  Real qc[1+NVAPOR];
  std::fill(qc, qc + 1 + NVAPOR, 0.);
  for (int n = 1 + NVAPOR; n <= NVAPOR; ++n)
    qc[1+(n-1)%NVAPOR] += prim[n] + 1.0E-10;  // prevent devide by 0

  Real lv = 0.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NVAPOR;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    lv += GetLatent(n,temp)*prim[ng]*ratio;
  }

  Real st = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) {
    int ng = 1 + (n-1)%NVAPOR;
    Real ratio = (prim[n] + 1.0E-10)/qc[ng];
    st += (prim[n] + prim[ng]*ratio)*(GetCpRatio(n) - 1.);
  }
  Real lv_ov_cpt = lv/(cpd*st*temp);

  Real chi = Rd_/cpd*qd/st;

  Real xv = 1.;
  for (int n = 1; n <= NVAPOR; ++n)
    xv += prim[n]/qd/mu_ratios_[n];

  Real pd = pres/xv;

  Real rh = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    Real eta = prim[n]/qd/mu_ratios_[n];
    Real pv = pres*eta/xv;
    int nc = n + NVAPOR;
    Real esat;
    if (n == AMMONIA_VAPOR_ID)
      esat = sat_vapor_p_NH3_BriggsS(temp);
    else if (n == WATER_VAPOR_ID)
      esat = sat_vapor_p_H2O_BriggsS(temp);
    else
      esat = SatVaporPresIdeal(temp/t3_[nc], p3_[nc], beta_[nc], delta_[nc]);
    rh *= pow(pv/esat, -eta*Rd_/(cpd*st));
  }

  return temp*pow(p0/pd, chi)*exp(lv_ov_cpt)*rh;
#else
  return PotentialTemp(prim, p0, pthermo);
#endif*/
}


#endif /* end of include guard THERMODYNAMICS_FUNCS_HPP */

