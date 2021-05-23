#ifndef THERMODYNAMICS_HPP
#define THERMODYNAMICS_HPP

// C/C++ headers
#include <cfloat>
#include <iosfwd>

// Athena headers
#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../athena.hpp"
#include "thermodynamic_funcs.hpp"
#include "vapors/ammonia_vapors.hpp"
#include "vapors/water_vapors.hpp"

class MeshBlock;
class ParameterInput;
template<typename T> class AthenaArray;

// dynamic variables are ordered in the array as the following
// 0: dry air (non-condensible)
// 1+iphase*NVAPOR:1+(iphase+1)*NVAPOR
// iphase = 0 - moist air (condensible)
// iphase = 1 - primary condensible species
// iphase = 2..(N-1) - all other condensible species

enum class Adiabat {reversible = 0, pseudo = 1, dry = 2, isothermal = 3};

class Thermodynamics {
  friend std::ostream& operator<<(std::ostream& os, Thermodynamics const& my);
public:
  // functions
  Thermodynamics(MeshBlock *pmb, ParameterInput *pin);
  ~Thermodynamics();

  //! ideal gas constant of dry air in J/(kg K)
  Real GetRd() const {return Rd_;}

  /*! ratio of specific heat capacity at constant volume
   * @param n index of vapor
   * @return $c_v$ of vapor n in units of J/(kg K) ratioed to that of dry air
   */
  Real GetCvRatio(int n) const {return rcv_[n];}

  /*! specific heat capacity at constant volume in J/(kg K)
   * @param n index of vapor
   */
  Real GetCv(int n) const {
    Real cvd = Rd_/(pmy_block_->peos->GetGamma() - 1.);
    return rcv_[n]*cvd;
  }

  /*! specific heat capacity at constant pressure in J/(kg K)
   * @param n index of vapor
   */
  Real GetCpRatio(int n) const {return rcp_[n];}

  /*! ratio of specific heat capacity at constant pressure
   * @param n index of vapor
   * @return $c_p$ of vapor n in units of J/(kg K) ratioed to that of dry air
   * */
  Real GetCp(int n) const {
    Real gamma = pmy_block_->peos->GetGamma();
    Real cpd = Rd_*gamma /(gamma - 1.);
    return rcp_[n]*cpd;
  }

  /*! temperature dependent latent heat of condensates in J/kg (positive)
   * @param n index of condensate
   * @param temp temperature (default to 0)
   */
  Real GetLatent(int n, Real temp = 0.) const {
    return latent_[n] - delta_[n]*Rd_/eps_[n]*temp;
  }

  //! ratio of molecular weights
  Real GetMassRatio(int n) const {return eps_[n];}

  //! beta parameter
  Real GetBeta(int n) const {return beta_[n];}

  //! delta parameter
  Real GetDelta(int n) const {return delta_[n];}

  //! saturation vapor pressure
  Real SatVaporPressure(Real temp, int n) const {
    if (n == AMMONIA_VAPOR_ID)
      return sat_vapor_p_NH3_BriggsS(temp);
    else if (n == WATER_VAPOR_ID)
      return sat_vapor_p_H2O_BriggsS(temp);
    else
      return SatVaporPresIdeal(temp/t3_[n], p3_[n], beta_[n], delta_[n]);
  }

  /*! Construct an adiabat with primary condensate
   * @param method = 0 - reversible adiabat \
   *               = 1 - pseudo adiabat \
   *               = 2 - an adiabat with no latent heat release
   */
  void ConstructAdiabat(Real **w, Real Ts, Real Ps,
    Real grav, Real dz, int len, Adiabat method, Real dTdz = 0.) const;

  // uhat is the molar interal energy defined as:
  // u^\hat = (q_d^\hat c_d^\hat T 
  //        + \sum_i q_i^\hat c_i^\hat T
  //        + \sum_{i,j} q_{ij}^\hat c_{ij}^\hat T
  //        - \sum_{i,j} q_{ij}^\hat \mu_{ij}^\hat)/R_d^\hat
  void UpdateTPConservingU(Real q[], Real rho, Real uhat) const;

  /*! adjust conserved variable to a sub-saturated state
   * @param u conserved variables
   */
  void SaturationAdjustment(AthenaArray<Real> &u) const;

  //void PolytropicIndex(AthenaArray<Real> &gamma, AthenaArray<Real> &w,
  //  int kl, int ku, int jl, int ju, int il, int iu) const;

  //! polytropic index $\gamma=c_p/c_v$
  template<typename T>
  Real GetPolytropicIndex(T w);

  //! Conserved variables to thermodynamic variables
  void ConservedToThermo(AthenaArray<Real> &q, AthenaArray<Real> const& u,
    int il, int iu, int jl, int ju, int kl, int ku) const;

  // Thermodynamic variables to conserved variables
  // density of dry air, momentum and total energy are not updated
  void ThermoToConserved(AthenaArray<Real> &u, AthenaArray<Real> const& q,
    int il, int iu, int jl, int ju, int kl, int ku) const;

  // template functions
  template<typename T>
  void Q2P(T w, Real const q[]) {
    cut_hat(w1_, q, eps_, Rd_);
    for (int n = 0; n < NHYDRO; ++n) w[n] = w1_[n];
  }

  template<typename T>
  void P2Q(Real q[], T const w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    put_hat(q, w1_, eps_, Rd_);
  }

  template<typename T>
  Real GetChi(T w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    Real gamma = pmy_block_->peos->GetGamma();
    Real tem[1] = {Temp(w)};
    update_gamma(gamma, rcp_, tem);
    Real qeps = q_eps(w1_, eps_);
    Real qsig = q_sig(w1_, rcp_);
    return (gamma - 1.)/gamma*qeps/qsig;
  }

  template<typename T>
  Real GetCp(T w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    Real gamma = pmy_block_->peos->GetGamma();
    Real tem[1] = {Temp(w)};
    update_gamma(gamma, rcp_, tem);
    Real qsig = q_sig(w1_, rcp_);
    return gamma/(gamma - 1.)*Rd_*qsig;
  }

  template<typename T>
  Real GetCv(T w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    Real gamma = pmy_block_->peos->GetGamma();
    Real tem[1] = {Temp(w)};
    update_gamma(gamma, rcp_, tem);
    Real qsig = q_sig(w1_, rcp_);
    return 1./(gamma - 1.)*Rd_*qsig;
  }

  template<typename T>
  Real Qeps(T w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    return q_eps(w1_, eps_);
  }

  template<typename T>
  Real GetTemp(T w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    return w[IPR]/(w[IDN]*Rd_*q_eps(w1_, eps_));
  }

  template<typename T>
  Real GetTheta(T w, Real p0) {
    Real chi = Chi(w);
    Real temp = Temp(w);
    return temp*pow(p0/w[IPR], chi);
  }

  template<typename T>
  Real GetThetaE(T prim, Real p0) {
#if (NVAPOR > 0)
    Real gamma = pmy_block_->peos->GetGamma();
    Real tem[1] = {Temp(prim)};
    update_gamma(gamma, const_cast<Real*>(rcp_), tem);
    Real cpd = Rd_*gamma/(gamma - 1.);
    Real temp = Temp(prim);
    Real pres = prim[IPR];

    Real qd = 1.;
    for (int n = 1; n < NMASS; ++n)
      qd -= prim[n];

    Real qc[1+NVAPOR];
    std::fill(qc, qc + 1 + NVAPOR, 0.);
    for (int n = 1 + NVAPOR; n < NMASS; ++n)
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
      xv += prim[n]/qd/eps_[n];

    Real pd = pres/xv;

    Real rh = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      Real eta = prim[n]/qd/eps_[n];
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
    return Theta(prim, p0);
#endif
  }

  template<typename T>
  Real GetMoistStaticEnergy(T prim, Real gz) {
    Real LE = 0.;
    for (int n = 1 + NVAPOR; n < NMASS; ++n)
      LE -= latent_[n]*prim[n];
    return Cp(prim)*Temp(prim) + LE + gz;
  }

  /*! Saturation surplus for vapors can be both positive and negative
   * positive value represents supersaturation
   * negative value represents saturation deficit
   */
  void CalculateSaturationSurplus(Real dq[], AthenaArray<Real> const& w) {
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
        // default to boilding (evaporate all)
        if (xs < 1.) {
          Real ys = xs/(1. - xs);
          dq[n] = std::max(dq[n], prim[n]*(1. - ys/yy)*sum/qtol);
        }
      }
    }
  }

  // Real last_time, dt;
private:
  MeshBlock *pmy_block_;
  Real ftol_;
  int max_iter_;
  Real w1_[NHYDRO];   /**< scratch array for storing primitive variables */

  Real Rd_;           /**< ideal gas constant of dry air */
  Real eps_[NMASS];   /**< ratio of mean molecular weights */
  Real rcp_[NMASS];   /**< ratio of heat capacity at constant pressure */
  Real beta_[NMASS];  /**< dimensionless latent heat */

  Real latent_[NMASS];  /**< latent heat in J/kg */
  Real delta_[NMASS];   /**< cp(gas) - cp(cloud) */
  Real rcv_[NMASS];     /**< ratio of heat capacity at constant volume */
  Real t3_[NMASS];      /**< triple point temperature */
  Real p3_[NMASS];      /**< triple point pressure */
};

#endif
