/** @file thermodynamics_impl.hpp
 * @brief Implementation of templated functions in class Thermodynamics
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday May 23, 2021 12:32:22 EDT
 * @bug No known bugs.
 */

#ifndef THERMODYNAMICS_IMPL_HPP
#define THERMODYNAMICS_IMPL_HPP

template<typename T>
Real Thermodynamics::GetPolytropicIndex(T w) {
  Real gamma = pmy_block_->peos->GetGamma();
  Real fsig = 1., feps = 1.;
  for (int n = 1; n <= NVAPOR; ++n) {
    fsig += w[n]*(rcv_[n] - 1.);
    feps += w[n]*(1./eps_[n] - 1.);
  }
  return 1. + (gamma - 1.)*feps/fsig;
}

template<typename T>
void Molar2Mass(T w, Real const q[]) {
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


#endif /* end of include guard THERMODYNAMICS_IMPL_HPP */

