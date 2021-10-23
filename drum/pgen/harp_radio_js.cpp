// C/C++ header
#include <ctime>
#include <iomanip>
#include <iostream>
#include <fstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/linalg.h"
#include "../thermodynamics/thermodynamic_funcs.hpp"
#include "../thermodynamics/molecules.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/microwave/mwr_absorbers.hpp"
#include "../inversion/inversion.hpp"
#include "../inversion/radio_data.hpp"

// molecules
enum {iH2O = 1, iNH3 = 2};
Real grav, P0, T0, Z0, Tmin, xHe, xCH4;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(4+NVAPOR);
  SetUserOutputVariableName(0, "temp", "temperature", "K");
  SetUserOutputVariableName(1, "theta", "potential temperature", "K");
  SetUserOutputVariableName(2, "thetav", "virtual potential temperature", "K");
  SetUserOutputVariableName(3, "mse", "moist static energy", "J/kg");
  for (int n = 1; n <= NVAPOR; ++n) {
    std::string name = "rh" + std::to_string(n);
    SetUserOutputVariableName(3+n, name.c_str(), "relative humidity");
  }
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;
  Real p1bar = 1.E5;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = PotentialTemp(phydro->w.at(k,j,i), p1bar, pthermo);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(phydro->w.at(k,j,i));
        // mse
        user_out_var(3,k,j,i) = MoistStaticEnergy(phydro->w, grav*pcoord->x1v(i),
          pthermo, ppart, k, j, i);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3+n,k,j,i) = RelativeHumidity(phydro->w.at(k,j,i), n, pthermo);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  Z0 = pin->GetOrAddReal("problem", "Z0", 0.);
  Tmin = pin->GetReal("problem", "Tmin");
  grav = - pin->GetReal("hydro", "grav_acc1");

  std::string planet = pin->GetOrAddString("job", "planet", "");
  Real latitude = pin->GetOrAddReal("job", "latitude", 0.);
  if (planet != "") // update gravity at specific latitude
    grav = GetGravity(planet.c_str(), latitude);
}

void RadiationBand::AddAbsorber(std::string name, std::string file, ParameterInput *pin)
{
  std::stringstream msg;

  xHe = pin->GetReal("problem", "xHe");
  xCH4 = pin->GetReal("problem", "xCH4");

  if (name == "mw_CIA") {
    pabs->AddAbsorber(MwrAbsorberCIA(this, xHe, xCH4));
  } else if (name == "mw_NH3") {
    pabs->AddAbsorber(MwrAbsorberNH3(this, {iNH3, iH2O}, xHe).SetModelBellottiSwitch());
    //pabs->AddAbsorber(MwrAbsorberNH3(this, {iNH3, iH2O}, xHe).SetModelHanley());
  } else if (name == "mw_H2O") {
    pabs->AddAbsorber(MwrAbsorberH2O(this, iH2O, xHe));
  } else {
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknow absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

// If you want to use real gas cp
/*void update_gamma(Real& gamma, Real const q[]) {
	Real T = q[IDN], cp_h2, cp_he, cp_ch4;
	if (T < 300.)
    cp_h2 = Hydrogen::cp_norm(T);
	else
    cp_h2 = Hydrogen::cp_nist(T);
  cp_he = Helium::cp_nist(T);
  cp_ch4 = Methane::cp_nist(T);

	Real cp_real = (1. - xHe - xCH4)*cp_h2 + xHe*cp_he + xCH4*cp_ch4;
  gamma = cp_real/(cp_real - Thermodynamics::Rgas);
}*/

// implemented in src/inversion/radio_observation_lnp.cpp
Real RadioObservationLnProb(Real *par, Real *val, int ndim, int nvalue, void *pdata);

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  ATHENA_LOG("harp_radio_js");
  std::stringstream msg;
  //ReadJunoMWRProfile("Juno_MWR_PJ1345689_m24-m16_avgcoeff.fits", coeff, cov);
  //ReadWriteGeminiTEXESMap("Gemini_TEXES_GRS_2017_product.dat", coeff, iNH3);

  // 1. construct a 1D pseudo-moist adiabat at (T0,P0,Z0)
  Real gamma = pin->GetReal("hydro", "gamma");

  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real x1rat = pmy_mesh->mesh_size.x1rat;

  Real dz, **w1, *z1, *p1, *t1;
  if (x1rat != 1.0) {
    dz = (x1max - x1min)*(x1rat - 1.)/(pow(x1rat, pmy_mesh->mesh_size.nx1) - 1.);
    dz = std::min(dz, dz*pow(x1rat, pmy_mesh->mesh_size.nx1))/2.;
  } else {
    dz = (x1max - x1min)/pmy_mesh->mesh_size.nx1/2.;
  }
  int nx1 = (int)((x1max - x1min)/dz);
  NewCArray(w1, nx1, NHYDRO+2*NVAPOR);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // 1.1 estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*(x1min - Z0);
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;

  for (int n = 1; n <= NVAPOR; ++n) {
    Real qv = pin->GetReal("problem", "qvapor" + std::to_string(n))/1.E3;
    w1[0][n] = qv;
  }
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  std::cout << "\t## Request T = " << T0 << " P = " << P0 << " at Z = " << Z0 << std::endl;
  while (iter++ < max_iter) {
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo);

    // 1.2 replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1; ++ii)
      if (pthermo->GetTemp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IPR] = w1[ii][IPR]*exp(-grav*(z1[i] - z1[ii])/(Rd*Tv));
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n <= NVAPOR; ++n)
        w1[i][n] = w1[ii][n];
    }

    // 1.3 find TP at z = Z0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->GetTemp(w1[i]);
    }
    p0 = interp1(Z0, p1, z1, nx1);
    t0 = interp1(Z0, t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0/p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0/p0 - 1.) < 1.E-4)) break;
    std::cout << "\t## Iteration #" << iter << ": " << "T = " << t0 << " P = " << p0 << std::endl;
  }

  if (iter > max_iter) {
    std::stringstream msg;
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0
        << std::endl << "P0 = " << p0;
    ATHENA_ERROR(msg);
  }

  // setup initial condition
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;

    // set gas concentration
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j)
          phydro->w(n,k,j,i) = buf[n];
  }

  // Apply boundary condition
  phydro->hbvar.SwapHydroQuantity(phydro->w, HydroBoundaryQuantity::prim);
  pbval->ApplyPhysicalBoundaries(0., 0.);

  // Calculate baseline radiation
  prad->CalculateRadiances(phydro->w, 0., ks, js, is, ie+1);

  //bool run_forward= pin->GetInteger("time", "nlim") > 0;
  bool run_forward= true;

  // Initialize objective function
  RadioData myradio = RadioData(this, pin, 
    grav, iH2O, iNH3, z1, p1, t1, nx1, run_forward);

  int nwalker = pin->GetInteger("problem", "nwalker");
  int ndim = 3*myradio.zfrac.size();  // location, T and NH3
  int nwave = prad->GetNumBands();

  // parameter array and output value array
  Real **par, **val;
  NewCArray(par, nwalker, ndim); 
  NewCArray(val, nwave, 3);

  int nsample = myradio.TpSample.size();
  for (int i = 0; i < nsample; ++i) {
    par[0][i] = myradio.zfrac[i];
    if (myradio.zfrac[i] < 0. || myradio.zfrac[i] > 1.) {
      msg << "### FATAL ERROR in ProblemGenerator" << std::endl 
          << "zfrac[" << i << "] should between [0.,1.]";
      ATHENA_ERROR(msg);
    }
    par[0][nsample+i] = myradio.TpSample[i];
    par[0][2*nsample+i] = myradio.NH3pSample[i];
  }

  /* Modify profile based on input
  //RadioObservationLnProb(*par, *val, ndim, nwave*3, &myradio);

  // copy modified model to baseline
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j < js+3; ++j)
        for (int i = is; i <= ie; ++i)
          phydro->w(n,k,j,i) = phydro->w(n,k,js+3,i);

  // recalculate baseline radiation
  prad->CalculateRadiances(phydro->w, 0., ks, js, is, ie + 1);

  // update baseline z,p,t
  for (int i = is; i <= ie; ++i) {
    myradio.z1[i-is] = pcoord->x1v(i);
    myradio.p1[i-is] = phydro->w(IPR,ks,js,i);
    myradio.t1[i-is] = pthermo->GetTemp(phydro->w.at(ks,js,i));
  }
  myradio.nx1 = ie - is + 1; // ie - is + 1 < nx1
  for (int i = 0; i < nsample; ++i)
    myradio.zdiv[i] = interp1(myradio.pdiv[i]*1.E5, myradio.z1, myradio.p1, myradio.nx1);

  // initialize random positions
  srand(time(NULL) + Globals::my_rank);

  for (int n = 0; n < nwalker; ++n) {
    for (int i = 0; i < ndim/3; ++i)
      par[n][i] = 1.*rand()/RAND_MAX;
    for (int i = ndim/3; i < 2*ndim/3; ++i)
      par[n][i] = (1.*rand()/RAND_MAX - 0.5)*myradio.Tstd;
    for (int i = 2*ndim/3; i < ndim; ++i)
      par[n][i] = (1.*rand()/RAND_MAX - 0.5)*myradio.NH3std;
  }

  pinvt->EnrollObjectives(RadioObservationLnProb, &myradio);
  pinvt->Initialize(par, nwalker, ndim, nwave*3);*/

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);

  FreeCArray(w1);
  FreeCArray(par);
  FreeCArray(val);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
