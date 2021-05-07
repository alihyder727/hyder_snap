// C/C++ header
#include <ctime>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../math/interpolation.h"
#include "../globals.hpp"
#include "../utils/utils.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../physics/physics.hpp"

// molecules
enum {iH2O = 1, iH2Oc = 2, iH2Op = 3};
Real grav, P0, T0, Tmin, prad, hrate, sponge_tau;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
  SetUserOutputVariableName(2, "thetav");
  SetUserOutputVariableName(3, "mse");
  SetUserOutputVariableName(4, "rh_h2o");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  Real dq[1+NVAPOR], rh;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->Temp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = pthermo->Theta(phydro->w.at(k,j,i), P0);
        // theta_v
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->Qeps(phydro->w.at(k,j,i));
        user_out_var(3,k,j,i) = pthermo->MSE(phydro->w.at(k,j,i), grav*pcoord->x1v(i));
        pthermo->SaturationSurplus(dq, phydro->w.at(k,j,i));
        rh = phydro->w(iH2O,k,j,i)/(phydro->w(iH2O,k,j,i) - dq[iH2O]);
        user_out_var(4,k,j,i) = rh;
      }
}

void Forcing(MeshBlock *pmb, Real const time, Real const dt,
    AthenaArray<Real> const &w, AthenaArray<Real> const &bcc, AthenaArray<Real> &u)
{
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;
  Thermodynamics *pthermo = pmb->pthermo;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real cv = pthermo->Cv(w.at(k,j,i));
        if (w(IPR,k,j,i) < prad) {
          u(IEN,k,j,i) += dt*hrate*w(IDN,k,j,i)*cv*(1. +
                          1.E-4*sin(2.*M_PI*rand()/RAND_MAX));
        }

        Real temp = pthermo->Temp(w.at(k,j,i));
        if (temp < Tmin) {
          u(IEN,k,j,i) += w(IDN,k,j,i)*cv*(Tmin - temp)/sponge_tau*dt;
        }
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  P0 = pin->GetReal("problem", "P0");
  T0 = pin->GetReal("problem", "T0");
  prad = pin->GetReal("problem", "prad");
  hrate = pin->GetReal("problem", "hrate")/86400.;  // K/day to K/s
  sponge_tau = pin->GetReal("problem", "sponge_tau");
  Tmin = pin->GetReal("hydro", "min_tem");

  EnrollUserExplicitSourceFunction(Forcing);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real qH2O = pin->GetReal("problem", "qH2O")/1.E3;
  Real gamma = pin->GetReal("hydro", "gamma");

  // construct a 1D pseudo-moist adiabat with given relative humidity
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
  NewCArray(w1, nx1, NHYDRO);
  z1 = new Real [nx1];
  p1 = new Real [nx1];
  t1 = new Real [nx1];

  // 1.1 estimate surface temperature and pressure
  Real Rd = pthermo->GetRd();
  Real cp = gamma/(gamma - 1.)*Rd;
  Real Ts = T0 - grav/cp*x1min;
  Real Ps = P0*pow(Ts/T0, cp/Rd);
  int max_iter = 200, iter = 0;

  w1[0][iH2O] = qH2O;
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  Real t0, p0;
  while (iter++ < max_iter) {
    pthermo->ConstructAdiabat(w1, Ts, Ps, grav, dz, nx1, Adiabat::pseudo);

    // 1.2 replace adiabatic atmosphere with isothermal atmosphere if temperature is too low
    int ii = 0;
    for (; ii < nx1; ++ii)
      if (pthermo->Temp(w1[ii]) < Tmin) break;
    Real Tv = w1[ii][IPR]/(w1[ii][IDN]*Rd);
    for (int i = ii; i < nx1; ++i) {
      w1[i][IPR] = w1[ii][IPR]*exp(-grav*(z1[i] - z1[ii])/(Rd*Tv));
      w1[i][IDN] = w1[i][IPR]/(Rd*Tv);
      for (int n = 1; n < NMASS; ++n)
        w1[i][n] = w1[ii][n];
    }

    // 1.3 find TP at z = 0
    for (int i = 0; i < nx1; ++i) {
      p1[i] = w1[i][IPR];
      t1[i] = pthermo->Temp(w1[i]);
    }
    p0 = interp1(0., p1, z1, nx1);
    t0 = interp1(0., t1, z1, nx1);

    Ts += T0 - t0;
    Ps *= P0/p0;
    if ((fabs(T0 - t0) < 0.01) && (fabs(P0/p0 - 1.) < 1.E-4)) break;
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
  int kl = block_size.nx3 == 1 ? ks : ks-NGHOST;
  int ku = block_size.nx3 == 1 ? ke : ke+NGHOST;
  int jl = block_size.nx2 == 1 ? js : js-NGHOST;
  int ju = block_size.nx2 == 1 ? je : je+NGHOST;
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO);
    buf[IVX] = buf[IVY] = buf[IVZ] = 0.;
    // set cloud to zero
    for (int n = 1+NVAPOR; n < NMASS; ++n)
      buf[n] = 0.;
    for (int n = 0; n < NHYDRO; ++n)
      for (int k = kl; k <= ku; ++k)
        for (int j = jl; j <= ju; ++j)
          phydro->w(n,k,j,i) = buf[n];
  }

  // setup open lower boundary
  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j) {
        for (int n = 0; n < NHYDRO; ++n)
          w1[0][n] = phydro->w(n,k,j,is);
        Real P1 = w1[0][IPR];
        Real T1 = pthermo->Temp(w1[0]);
        Real dz = pcoord->dx1f(is);
        // adiabatic extrapolate half a grid down
        pthermo->ConstructAdiabat(w1, T1, P1, grav, -dz/2., 2, Adiabat::reversible);
        for (int n = 0; n < NHYDRO; ++n)
          for (int i = 1; i <= NGHOST; ++i)
            phydro->w(n,k,j,is-i) = w1[1][n];
      }
  }

  srand(Globals::my_rank + time(0));
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  pphy->Initialize(phydro->w);

  FreeCArray(w1);
  delete[] z1;
  delete[] p1;
  delete[] t1;
}
