// C/C++ headers
#include <stdexcept>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../radiation/radiation.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/root.h"
#include "../math/linalg.h"
#include "inversion.hpp"
#include "gaussian_process.hpp"
#include "radio_data.hpp"

struct RootData {
  Real rho0;
  Real Tv0;
  Real dz;
  Real Rd;
  Real grav;
  Real chi;
  Real p0;
};

Real junomwr_root_func(Real Tv2, void *aux) {
  RootData *pd = static_cast<RootData*>(aux);
  Real& rho0 = pd->rho0;
  Real& Tv0 = pd->Tv0;
  Real& dz = pd->dz;
  Real& Rd = pd->Rd;
  Real& grav = pd->grav;
  Real& chi = pd->chi;
  Real& p0 = pd->p0;

  Real rho2 = rho0*(Tv0 - grav*dz/(2.*Rd))/(Tv2 + grav*dz/(2.*Rd));
  Real p2 = rho2*Rd*Tv2;
  if (rho2 > 0.)
    return Tv2 - Tv0*pow(p2/p0, chi);
  else
    return (Tv2 + Tv0)*log(Tv2/Tv0) + 2.*chi*grav*dz/Rd;
}

Real RadioObservationLnProb(Real *par, Real *val, int ndim, int nvalue, void *obj)
{
  // ######################   Load fitting data ############################# //
  RadioData *pobj = static_cast<RadioData*>(obj);
  std::vector<Real>& zfrac = pobj->zfrac;
  std::vector<Real>& TpSample = pobj->TpSample;
  std::vector<Real>& NH3pSample = pobj->NH3pSample;
  std::vector<Real> const& zdiv = pobj->zdiv;

  Real const& Tstd = pobj->Tstd;
  Real const& Tlen = pobj->Tlen;
  Real const& NH3std = pobj->NH3std;
  Real const& NH3len = pobj->NH3len;
  Real const& grav = pobj->grav;

  int const& iH2O = pobj->iH2O;
  int const& iNH3 = pobj->iNH3;

  int const& nx1 = pobj->nx1;
  int const& is = pobj->is;
  int const& ie = pobj->ie;
  int const& js = pobj->js;
  int const& je = pobj->je;
  int const& ks = pobj->ks;
  int const& ke = pobj->ke;

  Real const *z1 = pobj->z1;
  Real const *p1 = pobj->p1;
  Real const *t1 = pobj->t1;

  Real ***icov = pobj->icov;
  Real **target = pobj->target;

  std::vector<TPData> const& tpdata = pobj->tpdata;

  Thermodynamics *pthermo = pobj->pthermo;
  Hydro *phydro = pobj->phydro;
  Radiation *prad = pobj->prad;
  Coordinates *pcoord = pobj->pcoord;
  AthenaArray<Real> &w = phydro->w;

  // ######################   Forward modeling   ############################# //
  int nlayer = ie - is + 1;
  int nsample = TpSample.size();

  // 0..nsample-1 : zfrac, the first one is fixed
  // nsample..2*nsample-1 : TpSample
  // 2*nsample..3*nsample-1 : NH3pSample
  for (int i = 0; i < nsample; ++i) {
    zfrac[i] = par[i];
    if (zfrac[i] < 0. || zfrac[i] > 1.) return NAN;
    TpSample[i] = par[nsample + i];
    NH3pSample[i] = par[2*nsample + i];
  }

  // 1. copy baseline
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js+1; j <= je; ++j)
        for (int i = is; i <= ie; ++i)
          w(n,k,j,i) = w(n,k,js,i);

  // 2. locate height given pressure divides
  std::vector<Real> zlev(nsample), plev(nsample), tlev(nsample);

  if (pobj->IsTestRun())
    std::cout << "Sampling levels" << std::endl;
  for (int i = 0; i < nsample; ++i) {
    zlev[i] = zdiv[i]*(1.-zfrac[i]) + zdiv[i+1]*zfrac[i];
    plev[i] = interp1(zlev[i], p1, z1, nx1);
    tlev[i] = interp1(zlev[i], t1, z1, nx1);
    // report
    if (pobj->IsTestRun())
      printf("%12.4g km %12.4g bar %12.4g K\n", zlev[i]/1.E3, plev[i]/1.E5, tlev[i]);
  }

  // 3 calculate the covariance matrix of temperature
  Real *stdAll = new Real [nlayer];
  Real *stdSample = new Real [nsample];
  Real *Tp = new Real [nlayer];
  Real *NH3p = new Real [nlayer];

  // 4 calculate perturbed temperature profile
  for (int i = is; i <= ie; ++i)
    stdAll[i-is] = Tstd;
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = Tstd;

  gp_predict(SquaredExponential, Tp, pcoord->x1v.data() + is, stdAll, nlayer,
    TpSample.data(), zlev.data(), stdSample, nsample, Tlen);

  // fix boundary condition
  int ib = 0, it = nlayer - 1;
  while ((Tp[ib+1] < Tp[ib]) && (ib < nlayer)) ib++;
  for (int i = 0; i < ib; ++i) Tp[i] = Tp[ib]; 
  while ((Tp[it-1] > Tp[it]) && (it >= 0)) it--;
  for (int i = it+1; i < nlayer; ++i) Tp[i] = Tp[it]; 

  // 5 calculate perturbed ammonia profile
  for (int i = is; i <= ie; ++i)
    stdAll[i-is] = NH3std;
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = NH3std;
    
  gp_predict(SquaredExponential, NH3p, pcoord->x1v.data() + is, stdAll, nlayer,
    NH3pSample.data(), zlev.data(), stdSample, nsample, NH3len);

  // 6. save perturbed ammonia profile
  Real rho, Rd = pthermo->GetRd();
  int j1 = js+1, j2 = js+2, j3 = js+3;
  for (int i = is; i <= ie; ++i)
    w(iNH3,j1,i) += NH3p[i-is];

  // calculate new radiation if only ammonia is changed
  if (pobj->IsTestRun())
    prad->CalculateRadiances(w, 0., ks, j1, is, ie+1);

  // 7. save perturbed temperature profile
  for (int i = is; i <= ie; ++i) {
    Real temp = pthermo->GetTemp(w.at(js,i));
    w(IPR,j2,i) = w(IDN,js,i)*Rd*(temp + Tp[i-is])*pthermo->RovRd(w.at(js,i));
  }

  // calculate new radiation if only temperature is changed
  if (pobj->IsTestRun())
    prad->CalculateRadiances(w, 0., ks, j2, is, ie+1);

  // 8. rectify profile
  Real **w2;
  NewCArray(w2, 2, NHYDRO);
  std::fill(*w2, *w2 + 2*NHYDRO, 0.);
  w(IDN,j3,is) *= w(IPR,j3,is)/w(IPR,j2,is); 
  //std::cout << pthermo->Temp(w.at(js,is)) << " " << pthermo->Temp(w.at(j3,is)) << std::endl;
  for (int i = is; i < ie; ++i) {
    for (int n = 0; n < NHYDRO; ++n)
      w2[0][n] = w(n,j3,i);
    Real T0 = pthermo->GetTemp(w2[0]);
    Real P0 = w2[0][IPR];
    Real dz = pcoord->x1v(i+1) - pcoord->x1v(i);
    
    // propose an adiabatic move
    pthermo->ConstructAtmosphere(w2, T0, P0, grav, dz, 2, Adiabat::dry);
    Real T1 = pthermo->GetTemp(w2[1]);
    Real T2 = pthermo->GetTemp(w.at(js,i+1)) + Tp[i+1-is];
    bool isothermal_flag = false;

    if (T1 < 0.) {  // adiabat is not possible using isothermal instead
      w2[1][IDN] = w2[0][IDN];
      w2[1][IPR] = w2[0][IPR]*exp(-w2[0][IDN]*grav*dz/w2[0][IPR]);
      for (int n = 1; n <= NVAPOR; ++n)
        w2[1][n] = w2[0][n];
      T1 = pthermo->GetTemp(w2[1]);
      T2 = std::max(T1, T2);
      isothermal_flag = true;
    }

    Real Tv1 = w2[1][IPR]/(w2[1][IDN]*Rd);
    //std::cout << i << " " << T1 << " " << T2;
    if (T1 > T2) {  // potentially unstable
      //std::cout << " potentially unstable" << std::endl;
      Real eH2O = pthermo->GetMassRatio(iH2O);
      Real eNH3 = pthermo->GetMassRatio(iNH3);
      Real f = 1. + w2[1][iNH3]*(1./eNH3 - 1.);
      Real qmax = T1/T2*w2[0][iH2O] - eH2O/(eH2O - 1.)*(T1 - T2)/T2*f;
      Real Tv0 = w2[0][IPR]/(w2[0][IDN]*Rd);
      Real Tv2;
      //std::cout << Tv0 << " " << qmax << " " << w2[1][iH2O] << std::endl;
      if ((qmax > 0.) && (T2 > 0.)) { // has enough water
        w2[1][iH2O] = std::min(w2[1][iH2O], qmax);
        Tv2 = T2*pthermo->RovRd(w2[1]);
      } else {  // not enough water
        //std::cout << " not enough water" << std::endl;
        w2[1][iH2O] = 0.;
        //Tv2 = T2*pthermo->Qeps(w2[1]);
        RootData rd;
        rd.rho0 = w2[0][IDN];
        rd.Tv0 = Tv0;
        rd.dz = dz;
        rd.Rd = Rd;
        rd.grav = grav;
        rd.chi = pthermo->GetChi(w2[0]);
        rd.p0 = P0;
        //std::cout << "before Tv2=" << Tv2 << std::endl;
        if (Tv1 > grav*dz/Rd) {
          int err = root(Tv1/2, Tv0, 1.E-4, &Tv2, junomwr_root_func, &rd);
          if (err) {
            std::stringstream msg;
            msg << "### Root doesn't converge" << std::endl;
            msg << Tp[i+1-is] << " " << junomwr_root_func(Tv1/2, &rd) << " " 
                << junomwr_root_func(Tv0, &rd) << std::endl;
            std::cerr << msg.str() << std::endl;
            return NAN;
            //throw std::runtime_error(msg.str().c_str());
          }
        } else Tv2 = Tv0; // using isothermal
        // std::cout << "after Tv2=" << Tv2 << std::endl;
      }
      //std::cout << "before  "<< w2[1][IDN] << " " << w2[1][IPR] << " "
      //          << pthermo->Temp(w2[1]) << std::endl;
      w2[1][IDN] = w2[0][IDN]*(Tv0 - grav*dz/(2.*Rd))/(Tv2 + grav*dz/(2.*Rd));
      w2[1][IPR] = w2[1][IDN]*Rd*Tv2;
      //std::cout << "after " << w2[1][IDN] << " " << w2[1][IPR] << " "
      //          << pthermo->Temp(w2[1]) << std::endl;
    } else { // stable
      //std::cout << " (T2-T1)/dz=" << (T2-T1)/dz;
      //std::cout << " stable ";
      if (!isothermal_flag)
        pthermo->ConstructAtmosphere(w2, T0, P0, grav, dz, 2, Adiabat::dry, (T2-T1)/dz);
      //std::cout << pthermo->Temp(w2[1]) << " " << w2[1][IPR] << std::endl;
    }

    w2[1][iNH3] = std::max(0., w(iNH3,js,i+1) + NH3p[i+1-is]);
    Real dq[NHYDRO];
    pthermo->SaturationSurplus(dq, w2[1], VariableType::prim);
    if (dq[iNH3] > 0.) w2[1][iNH3] -= dq[iNH3];

    for (int n = 0; n < NHYDRO; ++n)
      w(n,j3,i+1) = w2[1][n];
  }

  // 9. calculate radiation
  prad->CalculateRadiances(w, 0., ks, j3, is, ie+1);

  // ######################   Probability   ############################# //
  Real **B, *b, *misfit, *rr;
  std::vector<Direction> out_dir = prad->GetOutgoingRays();
  int nangle = prad->GetNumBands();
  Real lnprob = 0.;

  NewCArray(B, nangle, 3);
  b = new Real [nangle];
  misfit = new Real [nangle];
  rr = new Real [nangle];
  
  // 10. solve coefficients using least square
  for (int i = 0; i < nangle; ++i) {
    B[i][0] = 1.;
    B[i][1] = 1. - out_dir[i].mu;
    B[i][2] = (1. - out_dir[i].mu)*(1. - out_dir[i].mu);
  }

  // 11. log likelihood
  RadiationBand *pband = prad->pband;
  int i = 0;
  while (pband != NULL) {
    for (int j = 0; j < nangle; ++j)
      b[j] = pband->btoa(j,ks,j3);
    leastsq(B, b, nangle, 3);

    for (int j = 0; j < 3; ++j) {
      misfit[j] = b[j] - target[i][j];
      val[i*3 + j] = b[j];
    }

    mvdot(rr, icov[i], misfit, 3, 3);
    lnprob += -0.5*vvdot(misfit, rr, 3);

    i++;
    pband = pband->next;
  }

  if (pobj->IsTestRun())
    std::cout << "lnprob1 = " << lnprob << std::endl;

  // 12. addition TP constraints
  Real *z2 = new Real [nlayer];
  Real *p2 = new Real [nlayer];
  Real *t2 = new Real [nlayer];
  for (int i = is; i <= ie; ++i) {
    z2[i-is] = pcoord->x1v(i);
    p2[i-is] = w(IPR,j3,i);
    t2[i-is] = pthermo->GetTemp(w.at(j3,i));
  }

  if (tpdata.size() > 0) {
    for (int i = 0; i < tpdata.size(); ++i) {
      Real tem = interp1(tpdata[i].P, t2, p2, nlayer);
      //std::cout << tpdata[i].P << " " << tpdata[i].T << " " << tem << std::endl;
      lnprob += -0.5*(tem-tpdata[i].T)*(tem-tpdata[i].T)/(tpdata[i].ERR*tpdata[i].ERR);
    }

    if (pobj->IsTestRun())
      std::cout << "lnprob2 = " << lnprob << std::endl;
  }

  // 13. log prior
  // 13.1 ammonia
  lnprob += gp_lnprior(SquaredExponential, NH3pSample.data(),
    zlev.data(), stdSample, nsample, NH3len);

  if (pobj->IsTestRun())
    std::cout << "lnprob3 = " << lnprob << std::endl;

  // 13.2 temperature, prior probability on proposed perturbation
  for (int i = 0; i < nsample; ++i)
    stdSample[i] = 100.*Tstd;
  lnprob += gp_lnprior(SquaredExponential, TpSample.data(),
    zlev.data(), stdSample, nsample, Tlen);

  // 13.3 prior probability on actual temperature at same pressure
  Real *told = new Real [nlayer];
  for (int i = 0; i < nsample; ++i) {
    plev[i] = interp1(zlev[i], p2, z2, nlayer);
    tlev[i] = interp1(zlev[i], t2, z2, nlayer);
    told[i] = interp1(plev[i], t1, p1, nx1);
    TpSample[i] = tlev[i] - told[i];
    stdSample[i] = Tstd;
  }
  lnprob += gp_lnprior(SquaredExponential, TpSample.data(),
    zlev.data(), stdSample, nsample, Tlen);

  if (pobj->IsTestRun())
    std::cout << "lnprob4 = " << lnprob << std::endl;

  delete[] stdAll;
  delete[] stdSample;
  delete[] Tp;
  delete[] NH3p;

  FreeCArray(B);
  delete[] b;
  delete[] misfit;
  delete[] rr;

  delete[] z2;
  delete[] p2;
  delete[] t2;
  delete[] told;
  FreeCArray(w2);

  return lnprob;
  // ######################     End     ############################# //
}
