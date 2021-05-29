/** @file thermodynamics.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Tuesday May 25, 2021 18:04:45 UTC
 * @bug No known bugs.
 */

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
  // members
  MeshBlock *pmy_block; /**< pointer to MeshBlock */

  // member functions
  Thermodynamics(MeshBlock *pmb, ParameterInput *pin);
  ~Thermodynamics() {}

  /*! Ideal gas constant of dry air in [J/(kg K)]
   * @return $R_d=\hat{R}/m_d$
   */
  Real GetRd() const {
    return Rd_;
  }

  /*! Ratio of specific heat capacity at constant volume
   * @param i index of the vapor
   * @return $c_{v,i}/c_{v,d}$
   */
  Real GetCvRatio(int i) const {
    return cv_ratios_[i];
  }

  /*! Specific heat capacity at constant volume
   *
   * $c_{v,d} = \frac{R_d}{\gamma_d - 1}$ \n
   * $c_{v,i} = \frac{c_{v,i}}{c_{v,d}}\times c_{v,d}$
   * @param i index of the vapor
   * @return $c_v$ [J/(kg K)]
   */
  Real GetCv(int i) const {
    Real cvd = Rd_/(pmy_block->peos->GetGamma() - 1.);
    return cv_ratios_[i]*cvd;
  }

  /*! Ratio of specific heat capacity at constant pressure
   * @param i index of the vapor
   * @return $c_{p,i}/c_{p,d}$
   */
  Real GetCpRatio(int i) const {
    return cp_ratios_[i];
  }

  /*! Specific heat capacity at constant pressure
   *
   * $c_{p,d} = \frac{\gamma_d}{\gamma_d - 1}R_d$ \n
   * $c_{p,i} = \frac{c_{p,i}}{c_{p,d}}\times c_{p,d}$
   * @param i index of the vapor
   * @return $c_p$ [J/(kg K)]
   * */
  Real GetCp(int i) const {
    Real gamma = pmy_block->peos->GetGamma();
    Real cpd = Rd_*gamma /(gamma - 1.);
    return cp_ratios_[i]*cpd;
  }

  /*! Temperature dependent specific latent heat of condensates at constant volume
   *
   * $L_{ij}(T) = L_{ij}^r - (c_{ij} - c_{p,i})\times(T - T^r)$
   * $= L_{ij}^r - \delta_{ij}R_i(T - T^r)$
   * @param j index of condensate
   * @param temp temperature (default to 0)
   * @return $L_{ij}(T)$ [J/kg]
   */
  Real GetLatent(int j, Real temp = 0.) const {
    return latent_[j] - delta_[j]*Rd_/mu_ratios_[j]*temp;
  }

  /*! Ratio of molecular weights
   * @param i index of the vapor or the condensate
   * @return $\epsilon_i=m_i/m_d$
   */
  Real GetMassRatio(int i) const {
    return mu_ratios_[i];
  }

  /*! $\beta$ parameter of a condensate
   *
   * $\beta_{ij} = \frac{\Delta U_{ij}}{R_i T^r}$ \n
   * $\Delta U_{ij}$ is the difference in internal energy between the vapor $i$ 
   * and the condensate $j$ \n
   * $R_i=\hat{R}/m_i=R_d/\epsilon_i$ \n
   * $T^r$ is the triple point temperature
   * @param j index of the condensate
   * @return $\beta_{ij}$
   */
  Real GetBeta(int j) const {
    return beta_[j];
  }

  /*! $\delta$ parameter of a condensate
   *
   * $\delta_{ij} = \frac{\Delta c_{ij}}{R_i}$ \n
   * $\Delta c_{ij} = c_{ij} - c_{p,j}$ is the difference in specific heat capacity 
   * at constant temperature.\n
   * $c_{ij}$ is the heat capacity of condensate $j$ of vapor $i$
   * @param j index of the condensate
   * @return $\delta_{ij}$
   */
  Real GetDelta(int j) const {
    return delta_[j];
  }

  /*! Construct an 1d atmosphere
   * @param method = 0 - reversible adiabat \n
   *               = 1 - pseudo adiabat \n
   *               = 2 - an adiabat with no latent heat release
   */
  void ConstructAtmosphere(Real **w, Real Ts, Real Ps,
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
  //void SaturationAdjustment(AthenaArray<Real> &u, AthenaArray<Real> &c) const;

  /*! Conserved variables to thermodynamic variables
   * @param q thermodynamic variables
   * @param u conserved variables
   * @param il start index in dim 1
   * @param iu end index in dim 1
   * @param jl start index in dim 2
   * @param ju end index in dim 2
   * @param kl start index in dim 3
   * @param ku end index in dim 3
   */
  void ConservedToThermo(AthenaArray<Real> &q, AthenaArray<Real> const& u,
    int il, int iu, int jl, int ju, int kl, int ku) const;

  // Thermodynamic variables to conserved variables
  // density of dry air, momentum and total energy are not updated
  void ThermoToConserved(AthenaArray<Real> &u, AthenaArray<Real> const& q,
    int il, int iu, int jl, int ju, int kl, int ku) const;

  //void PolytropicIndex(AthenaArray<Real> &gamma, AthenaArray<Real> &w,
  //  int kl, int ku, int jl, int ju, int il, int iu) const;

  //! polytropic index $\gamma=c_p/c_v$
  template<typename T>
  Real GetGamma(T w) {
    Real gamma = pmy_block->peos->GetGamma();
    Real fsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      fsig += w[n]*(cv_ratios_[n] - 1.);
      feps += w[n]*(1./mu_ratios_[n] - 1.);
    }
    return 1. + (gamma - 1.)*feps/fsig;
  }

  // template functions
  template<typename T>
  void ChemicalToPrimitive(T w, Real const q[]) {
    molar_to_mass(w1_, q, mu_ratios_, Rd_);
    for (int n = 0; n < NHYDRO; ++n) w[n] = w1_[n];
  }

  template<typename T>
  void PrimitiveToChemical(Real q[], T const w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    mass_to_molar(q, w1_, mu_ratios_, Rd_);
  }

  template<typename T>
  Real RovRd(T w) {
    Real feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n)
      feps += w[n]*(1./mu_ratios_[n] - 1.);
    return feps;
  }

  //! Temperature
  template<typename T>
  Real GetTemp(T w) {
    return w[IPR]/(w[IDN]*Rd_*RovRd(w));
  }

  template<typename T>
  Real GetPres(T u) {
    Real gm1 = pmy_block->peos->GetGamma() - 1;
    Real fsig = 0., feps = 0., rho = u[IDN];
    for (int n = 0; n <= NVAPOR; ++n) {
      rho += u[n];
      fsig += u[n]*cv_ratios_[n];
      feps += u[n]/mu_ratios_[n];
    }
    Real KE = 0.5*(u[IM1]*u[IM1] + u[IM2]*u[IM2] + u[IM3]*u[IM3])/rho;
    return gm1*(u[IEN] - KE)*feps/fsig;
  }

  template<typename T>
  Real GetChi(T w) {
    for (int n = 0; n < NHYDRO; ++n) w1_[n] = w[n];
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1., feps = 1.;
    for (int n = 1; n <= NVAPOR; ++n) {
      qsig += w[n]*(cp_ratios_[n] - 1.);
      feps += w[n]*(1./mu_ratios_[n] - 1.);
    }
    return (gamma - 1.)/gamma*feps/qsig;
  }

  template<typename T>
  Real GetCp(T w) {
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1.;
    for (int n = 1; n <= NVAPOR; ++n)
      qsig += w[n]*(cp_ratios_[n] - 1.);
    return gamma/(gamma - 1.)*Rd_*qsig;
  }

  template<typename T>
  Real GetCv(T w) {
    Real gamma = pmy_block->peos->GetGamma();
    Real tem[1] = {GetTemp(w)};
    update_gamma(gamma, tem);
    Real qsig = 1.;
    for (int n = 1; n <= NVAPOR; ++n)
      qsig += w[n]*(cv_ratios_[n] - 1.);
    return 1./(gamma - 1.)*Rd_*qsig;
  }

  //! Potential temperature
  template<typename T>
  Real GetTheta(T w, Real p0) {
    Real chi = GetChi(w);
    Real temp = GetTemp(w);
    return temp*pow(p0/w[IPR], chi);
  }

  /*! Saturation surplus for vapors can be both positive and negative
   * positive value represents supersaturation \n
   * negative value represents saturation deficit
   */
  template<typename T>
  void SaturationSurplus(Real dw[], T v, VariableType vtype = VariableType::prim) {
    // mass to molar mixing ratio
    if (vtype == VariableType::prim) {
      Real sum = 1.;
      for (int n = 1; n <= NVAPOR; ++n) {
        w1_[n] = v[n]/mu_ratios_[n];
        sum += v[n]*(1./mu_ratios_[n] - 1.);
      }
      for (int n = 1; n <= NVAPOR; ++n)
        w1_[n] /= sum;

      w1_[IDN] = GetTemp(v);
      w1_[IPR] = v[IPR];
    } else if (vtype == VariableType::cons) {
      Real sum = 0., feps = 0.;
      for (int n = 0; n <= NVAPOR; ++n) {
        w1_[n] = v[n]/mu_ratios_[n];
        sum += w1_[n];
      }
      for (int n = 0; n <= NVAPOR; ++n)
        w1_[n] /= sum;
      w1_[IPR] = GetPres(v);
      w1_[IDN] = w1_[IPR]/(sum*Rd_);
    } else {  // vtype == VariableType::chem
      for (int n = 0; n < NHYDRO+2*NVAPOR; ++n)
        w1_[n] = v[n];
    }

    for (int iv = 1; iv <= NVAPOR; ++iv) { 
      int nc = w1_[IDN] > t3_[iv] ? iv + NVAPOR : iv + 2*NVAPOR;
      int ic = NHYDRO - NVAPOR + nc - 1;
      Real rate = VaporCloudEquilibrium(w1_, iv, ic, t3_[iv], p3_[iv], 
          0., beta_[nc], delta_[nc], true);
      dw[iv] = rate/w1_[iv]*v[iv];
    }
  }

  //! Relative humidity
  template<typename T>
  Real GetRelativeHumidity(T w, int iv) {
    SaturationSurplus(dw_, w);
    return w[iv]/(w[iv] - dw_[iv]);
  }

  //! Equivalent potential temperature
  //template<typename T>
  //Real GetThetaE(T prim, Real p0);

private:
  Real ftol_;
  int max_iter_;

  //! scratch array for storing variables
  Real w1_[NHYDRO+2*NVAPOR];

  //! scratch array for storing variables
  Real dw_[1+NVAPOR];

  // read from inputs
  //! ideal gas constant of dry air in J/kg
  Real Rd_;
  //! ratio of mean molecular weights
  Real mu_ratios_[1+3*NVAPOR];
  //! ratio of specific heat capacities at constant pressure
  Real cp_ratios_[1+3*NVAPOR];
  /*! dimensionless latent heat
   *$\beta_{ij} == \frac{L_{ij}^r + \Delta c_{ij}T^r}{R_i T^r} $
   */
  Real beta_[1+3*NVAPOR];
  //! triple point temperature [K]
  Real t3_[1+NVAPOR];
  //! triple point pressure [pa]
  Real p3_[1+NVAPOR];

  // calculated quantities
  /*! latent heat in J/kg
   *$L_i^r+\Delta c_{ij}T^r == \beta_{ij}\frac{R_d}{\epsilon_i}T^r$
   */
  Real latent_[1+3*NVAPOR];

  /*! dimensionless differences in specific heat capacity
   *$(c_{ij} - c_{p,i})/R_i$
   */
  Real delta_[1+3*NVAPOR];

  //! ratio of specific heat capacities at constant volume
  Real cv_ratios_[1+3*NVAPOR];
};

//#if (NVAPOR > 0)
//  #include "get_theta_e.hpp"
//#endif

#endif
