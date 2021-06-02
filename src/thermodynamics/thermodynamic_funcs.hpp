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
  Real cpd = pthermo->GetCp(IDN);
  Real rho = w[IDN];
  for (int n = 0; n < 2*NVAPOR; ++n) {
    rho += c[n];
    IE += c[n]*(pthermo->GetCp(n)*temp - pthermo->GetLatent(1+NVAPOR+n));
  }
  return IE/rho + gz;
}

//! Relative humidity
template<typename T>
Real RelativeHumidity(T w, int iv, Thermodynamics *pthermo) {
  Real dw_[1+NVAPOR];
  pthermo->SaturationSurplus(dw_, w);
  return w[iv]/(w[iv] - dw_[iv]);
}

//! Equivalent potential temperature
//template<typename T>
//Real EquivalentPotentialTemp(T prim, Real p0, Thermodynamics *pthermo);


#endif /* end of include guard THERMODYNAMICS_FUNCS_HPP */

