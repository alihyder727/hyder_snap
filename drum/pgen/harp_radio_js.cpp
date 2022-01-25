/* -------------------------------------------------------------------------------------
 * SNAP Example Program
 *
 * Contributer:
 * Cheng Li, University of Michigan
 *
 * Year: 2022
 * Contact: chengcli@umich.edu
 * Reference: TBD
 * -------------------------------------------------------------------------------------
 */

// @sect3{Include files}

// C/C++ header
#include <iostream>

// Athena++ headers
// These input files are the same as those in the @ref straka. No additional comments
// are needed.
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamic_funcs.hpp"

// Here we include more input files
#include "../thermodynamics/molecules.hpp"
#include "../globals.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/linalg.h"
#include "../debugger/debugger.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/microwave/mwr_absorbers.hpp"
#include "../inversion/inversion.hpp"
#include "../inversion/radio_observation.hpp"

// @sect3{Preamble}

// molecules
enum {iH2O = 1, iNH3 = 2};
Real grav, P0, T0, Tmin, xHe, xCH4, rdlnTdlnP;

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

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = PotentialTemp(phydro->w.at(k,j,i), P0, pthermo);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(phydro->w.at(k,j,i));
        // mse
        user_out_var(3,k,j,i) = MoistStaticEnergy(phydro->w.at(k,j,i), grav*pcoord->x1v(i), pthermo);
        for (int n = 1; n <= NVAPOR; ++n)
          user_out_var(3+n,k,j,i) = RelativeHumidity(phydro->w.at(k,j,i), n, pthermo);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("hydro", "reference_pressure");
  T0 = pin->GetReal("problem", "T0");
  Tmin = pin->GetReal("problem", "Tmin");
  rdlnTdlnP = pin->GetOrAddReal("problem", "rdlnTdlnP", 1.);

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
    pabs->AddAbsorber(MwrAbsorberNH3(this, {iNH3, iH2O}, xHe).SetModelHanley());
  } else if (name == "mw_H2O") {
    pabs->AddAbsorber(MwrAbsorberH2O(this, iH2O, xHe));
  } else {
    msg << "### FATAL ERROR in RadiationBand::AddAbsorber"
        << std::endl << "unknow absorber: '" << name <<"' ";
    ATHENA_ERROR(msg);
  }
}

#ifdef REAL_GAS_CP
void update_gamma(Real& gamma, Real const q[]) {
  Real T = q[IDN], cp_h2, cp_he, cp_ch4;
  if (T < 300.)
    cp_h2 = Hydrogen::cp_norm(T);
  else
    cp_h2 = Hydrogen::cp_nist(T);
  cp_he = Helium::cp_nist(T);
  cp_ch4 = Methane::cp_nist(T);

  Real cp_real = (1. - xHe - xCH4)*cp_h2 + xHe*cp_he + xCH4*cp_ch4;
  gamma = cp_real/(cp_real - Thermodynamics::Rgas);
}
#endif

// @sect3{Initial condition}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  static_assert(HYDROSTATIC, "This problem requires turning on hydrostatic option");
  pdebug->Enter("ProblemGenerator: harp_radio_js");
  std::stringstream &msg = pdebug->msg;
  //ReadJunoMWRProfile("Juno_MWR_PJ1345689_m24-m16_avgcoeff.fits", coeff, cov);
  //ReadWriteGeminiTEXESMap("Gemini_TEXES_GRS_2017_product.dat", coeff, iNH3);

  // Construct a 1D pseudo-moist adiabat at (T0,P0)
  Real gamma = pin->GetReal("hydro", "gamma");
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;

  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  Real H0 = phydro->scale_height;

  // estimate surface temperature and pressure
  // coordinate range needs to be adjusted for hydrostatic model
  Real Ps = P0*exp(-x1min/H0);
  Real Ts = T0*pow(Ps/P0,Rd/cp);
  Real dlnp, **w1, *z1, *t1;
  dlnp = (x1max - x1min)/pmy_mesh->mesh_size.nx1/(2.*H0);
  int nx1 = (int)((x1max - x1min)/(H0*dlnp));
  NewCArray(w1, nx1, NHYDRO+2*NVAPOR);
  z1 = new Real [nx1];
  t1 = new Real [nx1];

  int max_iter = 200, iter = 0;

  for (int n = 1; n <= NVAPOR; ++n) {
    Real qv = pin->GetReal("problem", "qvapor" + std::to_string(n))/1.E3;
    w1[0][n] = qv;
  }
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + H0*dlnp;

  Real t0;
  if (Globals::my_rank == 0)
    std::cout << "- request T = " << T0 << " at P = " << P0 << std::endl;
  while (iter++ < max_iter) {
    pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, -dlnp, nx1, Adiabat::pseudo, rdlnTdlnP);

    // Replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1; ++ii)
      if (pthermo->GetTemp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n <= NVAPOR; ++n)
        w1[i][n] = w1[ii][n];
    }

    // Find T at p = p0
    for (int i = 0; i < nx1; ++i) {
      t1[i] = pthermo->GetTemp(w1[i]);
    }
    t0 = interp1(0, t1, z1, nx1);

    Ts += T0 - t0;
    if (Globals::my_rank == 0)
      std::cout << "- iteration #" << iter << ": " << "T = " << t0  << std::endl;
    if (fabs(T0 - t0) < 0.01) break;
  }

  if (iter > max_iter) {
    msg << "### FATAL ERROR in problem generator"
        << std::endl << "maximum iteration reached."
        << std::endl << "T0 = " << t0;
    ATHENA_ERROR(msg);
  }

  // Change to log quantity
  for (int i = 0; i < nx1; ++i)
    for (int n = 0; n < NHYDRO+2*NVAPOR; ++n) {
      if (n == IVX || n == IVY || n == IVZ)
        w1[i][n] = 0.;
      else
        w1[i][n] = log(w1[i][n]);
    }

  // setup initial condition
  // Z = -H0*log(P/P0)
  // log(P) = log(P0) - Z/H0
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    // set gas concentration and velocity
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = ks; k <= ke; ++k)
        for (int j = js; j <= je; ++j) {
          if (n == IVX || n == IVY || n == IVZ)
            phydro->w(n,k,j,i) = 0.;
          else
            phydro->w(n,k,j,i) = exp(buf[n]);
        }

  }

  // Apply boundary condition
  phydro->hbvar.SwapHydroQuantity(phydro->w, HydroBoundaryQuantity::prim);
  pbval->ApplyPhysicalBoundaries(0., 0.);

  // read profile updates from input
  std::vector<Real> TpSample, XpSample;
  std::vector<int> ix = {0, iNH3};

  Real Tstd = pin->GetReal("inversion", "Tstd");
  Real Tlen = pin->GetReal("inversion", "Tlen")*1.E3; // km -> m
  Real Xstd = pin->GetReal("inversion", "Xstd")*1.E-3;  // g/kg -> kg/kg
  Real Xlen = pin->GetReal("inversion", "Xlen")*1.E3; // km -> m
  TpSample = Vectorize<Real>(pin->GetString("problem", "Tp").c_str());
  XpSample = Vectorize<Real>(pin->GetString("problem", "NH3p").c_str());
  int nsample = pinvt->pradio->plevel.size();

  // add inversion boundary
  TpSample.insert(TpSample.begin(), 0.);
  TpSample.push_back(0.);
  XpSample.insert(XpSample.begin(), 0.);
  XpSample.push_back(0.);

  // update reference model
  for (int i = 0; i < nsample; ++i)
    XpSample[i] *= 1.E-3; // g/kg -> kg/kg

  // Revise atmospheric profiles and store them in position (ks,je).
  update_atm_profiles(this, ks, pinvt->pradio->plevel.data(), TpSample.data(), XpSample.data(),
      nsample, ix, Tstd, Tlen, Xstd, Xlen);

  // Copy revised profile to baseline position (ks,je) -> (ks,js).
  for (int n = 0; n < NHYDRO; ++n)
    for (int i = is; i <= ie; ++i)
      phydro->w(n,ks,js,i) = phydro->w(n,ks,je,i);

  // Calculate microwave radiation at (ks,js) at time t = 0
  // This value maintains for all future time steps (t > 0,ks,js)
  msg << "- running initial RT model 0" << std::endl;
  prad->CalculateRadiances(phydro->w, 0., ks, js, is, ie+1);

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  // Leave debug stack and print out all debug info.
  pdebug->Leave();

  // Clean up temporary memory
  FreeCArray(w1);
  delete[] z1;
  delete[] t1;
}
