// C/C++ headers
#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <regex>

// Athena++ headers
#include "../utils/utils.hpp"
#include "thermodynamic_funcs.hpp"
#include "vapors/ammonia_vapors.hpp"
#include "vapors/water_vapors.hpp"

Real dlnTdlnP(Real const q[], int const isat[], 
  Real const rcp[], Real const beta[], Real const delta[], Real const t3[], Real gamma)
{
  Real cp_hat = gamma/(gamma - 1.)*q_sig(q, rcp)/q_gas(q);

  Real xd = 1.;
  for (int n = 1; n < NMASS; ++n)
    xd -= q[n];

  Real c1 = 0., c2 = 0., c3 = 0.;
  for (int n = 1; n <= NVAPOR; ++n) {
    if (isat[n] > 0) {
      int nc = n + NVAPOR;
      Real latent = beta[nc]*t3[nc]/q[IDN] - delta[nc];
      c1 += q[n]/xd*latent;
      c2 += q[n]/xd*latent*latent;
    }
    c3 += q[n]/xd;
  }

  return (1. + c1)/(cp_hat + (c2 + c1*c1)/(1. + c3));
}

void __attribute__((weak)) update_gamma(Real& gamma, Real rcp[], Real const q[])
{}

// alpha = L/cv evaluated at current temperature
Real GasCloudIdeal(Real const q[], int iv, int ic,
  Real t3, Real p3, Real alpha, Real beta, Real delta)
{
  Real xv = q[iv];
  Real xc = q[ic];
  Real t = q[IDN]/t3;
  Real s;
  if (iv == AMMONIA_VAPOR_ID)
    s = sat_vapor_p_NH3_BriggsS(q[IDN]);
  else if (iv == WATER_VAPOR_ID)
    s = sat_vapor_p_H2O_BriggsS(q[IDN]);
  else
    s = SatVaporPresIdeal(t,p3,beta,delta);
  s /= q[IPR];

  // if saturation vapor pressure is larger than the total pressure
  // evaporate all condensates
  if (s > 1.) return -xc;

  Real g = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n) g -= q[n];
  g -= xv;

  Real s1 = s/(1. - s);
  Real rate = (xv - s1*g)/(1. + alpha*g*(beta/t - delta)*s1/(1. - s));

  // condensate at most xv vapor
  if (rate > 0.) rate = std::min(rate, xv);

  // evaporate at most xc cloud
  if (rate < 0.) rate = std::min(0., std::max(rate, -xc));
  return rate;
}

Real q_gas(Real const w[])
{
  Real fgas = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    fgas -= w[n];
  return fgas;
}

Real q_eps(Real const w[], Real const eps[])
{
  Real feps = 1.;
  for (int n = 1 + NVAPOR; n < NMASS; ++n)
    feps -= w[n];
  for (int n = 1; n <= NVAPOR; ++n)
    feps += w[n]*(1./eps[n] - 1.);
  return feps;
}

Real q_sig(Real const w[], Real const sig[])
{
  Real fsig = 1.;
  for (int n = 1; n < NMASS; ++n)
    fsig += w[n]*(sig[n] - 1.);
  return fsig;
}

Real qhat_eps(Real const q[], Real const eps[])
{
  Real feps = 1.;
  for (int n = 1; n < NMASS; ++n)
    feps += q[n]*(eps[n] - 1.);
  return q_gas(q)/feps;
}

void mass2molar(Real q[], Real const w[], Real const eps[], Real Rd)
{
  // set molar mixing ratio
  Real sum = 1.;
  for (int n = 1; n < NMASS; ++n) {
    q[n] = w[n]/eps[n];
    sum += w[n]*(1./eps[n] - 1.);
  }
  for (int n = 1; n < NMASS; ++n)
    q[n] /= sum;

  // set pressure and temperature
  q[IPR] = w[IPR];
  q[IDN] = w[IPR]/(w[IDN]*Rd*q_eps(w, eps));
}

void molar2mass(Real w[], Real const q[], Real const eps[], Real Rd)
{
  // set mass mixing ratio
  Real sum = 1.;
  for (int n = 1; n < NMASS; ++n) {
    w[n] = q[n]*eps[n];
    sum += q[n]*(eps[n] - 1.); 
  }
  for (int n = 1; n < NMASS; ++n)
    w[n] /= sum;

  // set pressure and density
  w[IPR] = q[IPR];
  w[IDN] = q[IPR]/(q[IDN]*Rd*q_eps(w, eps));
}

bool read_idealized_parameters(std::string fname, Real &pmin, Real &pmax, AthenaArray<Real> &dTdz,
  Real* &pmin_q, Real* &pmax_q, AthenaArray<Real>* &dqdz, std::vector<int> &mindex) 
{
  std::stringstream msg;
  if (!FileExists(fname)) {
    msg << "### FATAL ERROR in read_idealized_parameters. File " << fname 
        << " doesn't exist.";
    throw std::runtime_error(msg.str().c_str());
  }

  std::ifstream inp(fname.c_str(), std::ios::in);
  bool allocated = false;
  char c = '#';
  std::string line, ss = "";
  std::vector<std::string> arr;

  while (!inp.eof() && c == '#') {
    std::getline(inp, line);
    // removing leading space
    line = std::regex_replace(line, std::regex("^ +"), ss);
    c = line[0];
  }

  while (!inp.eof() && line.empty()) std::getline(inp, line);

  if (!inp.eof()) {
    arr = Vectorize<std::string>(line.c_str());
    int nq = arr.size() - 1;
    if (nq > 0) {
      mindex.resize(nq);
      for (int i = 0; i < nq; ++i)
        mindex[i] = std::stoi(arr[i+1]); // copy molecule index
      pmin_q = new Real[nq];
      pmax_q = new Real[nq];
      dqdz = new AthenaArray<Real>[nq];
      allocated = true;
    }

    // Read temperature parameters
    int nrows, ncols;
    inp >> pmin >> pmax >> nrows >> ncols;
    //std::cout << pmin << " " << pmax << " " << nrows << " " << ncols << std::endl;
    dTdz.NewAthenaArray(nrows, ncols);
    for (int i = 0; i < nrows; ++i)
      for (int j = 0; j < ncols; ++j)
        inp >> dTdz(i,j);

    // Read mixing ratios 
    for (int n = 0; n < nq; ++n) {
      inp >> pmin_q[n] >> pmax_q[n] >> nrows >> ncols;
      dqdz[n].NewAthenaArray(nrows, ncols);
      for (int i = 0; i < nrows; ++i)
        for (int j = 0; j < ncols; ++j)
          inp >> dqdz[n](i,j);
    }
  }

  return allocated;
}
