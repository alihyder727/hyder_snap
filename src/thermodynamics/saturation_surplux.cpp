/** @file saturation_surplux.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Saturday May 22, 2021 23:26:03 EDT
 * @bug No known bugs.
 */

#include "thermodynamics.hpp"

void CalculateSaturationSurplus(Real dq[], T prim) {
  Real temp = Temp(prim);

  // mass to molar mixing ratio
  Real sum = 1.;
  for (int n = 1; n < NMASS; ++n)
    sum -= prim[n]; // right now sum = qd
  dq[0] = sum;
  for (int n = 1; n <= NVAPOR; ++n) {
    dq[n] = prim[n]/eps_[n];
    sum += dq[n];
  }
  Real qtol = sum;
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    qtol += prim[n]/eps_[n];

  for (int n = 0; n <= NVAPOR; ++n)
    dq[n] /= sum;
  // Saturation surplus for vapors can be both positive and negative
  // positive value represents supersaturation
  // negative value represents saturation deficit

  for (int n = 1; n <= NVAPOR; ++n) {
    Real q = dq[n];
    Real yy = q/(1. - q);
    dq[n] = -1.0E10;
    for (int m = 1; m < NPHASE; ++m) {
      int nc = n + m*NVAPOR;
      Real svp;
      if (n == AMMONIA_VAPOR_ID)
        svp = sat_vapor_p_NH3_BriggsS(temp);
      else if (n == WATER_VAPOR_ID)
        svp = sat_vapor_p_H2O_BriggsS(temp);
      else
        svp = SatVaporPresIdeal(temp/t3_[nc], p3_[nc], beta_[nc], delta_[nc]);
      Real xs = svp/prim[IPR];

      // default to boiling (evaporate all)
      if (xs < 1.) {
        Real ys = xs/(1. - xs);
        dq[n] = std::max(dq[n], prim[n]*(1. - ys/yy)*sum/qtol);
      }
    }
  }
}

