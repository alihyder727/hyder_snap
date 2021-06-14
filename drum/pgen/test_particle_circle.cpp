//! \file test_thermodynamics.cpp
//  \brief Problem generator for testing thermodynamics

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../particles/particles.hpp"
#include "../utils/utils.hpp"
#include "../math/interpolation.h"
#include "../math/core.h"
#include "../thermodynamics/thermodynamic_funcs.hpp"

enum {iH2O = 1, iH2Oc = 2, iH2Op = 3};

Real p0, grav, omega;

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(5);
  SetUserOutputVariableName(0, "temp", "temperature", "K");
  SetUserOutputVariableName(1, "theta", "potential temperature", "K");
  SetUserOutputVariableName(2, "theta_v", "virtual potential temperature", "K");
  SetUserOutputVariableName(3, "mse", "moist static energy", "J/kg");
  SetUserOutputVariableName(4, "rh", "relative humidity");
  //SetUserOutputVariableName(4, "theta_e", "equivalent potential temperature", "K");
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        user_out_var(0,k,j,i) = pthermo->GetTemp(phydro->w.at(k,j,i));
        user_out_var(1,k,j,i) = PotentialTemp(phydro->w.at(k,j,i), p0, pthermo);
        user_out_var(2,k,j,i) = user_out_var(1,k,j,i)*pthermo->RovRd(phydro->w.at(k,j,i));
        user_out_var(3,k,j,i) = MoistStaticEnergy(phydro->w.at(k,j,i), ppart->c.at(k,j,i),
          grav*pcoord->x1v(i), pthermo);
        user_out_var(4,k,j,i) = RelativeHumidity(phydro->w.at(k,j,i), iH2O, pthermo);
        //user_out_var(4,i) = pthermo->GetThetaE(phydro->w.at(i), p0);
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  grav = - pin->GetReal("hydro", "grav_acc1");
  p0 = pin->GetReal("problem", "p0");
  omega = pin->GetOrAddReal("problem", "omega", 0.);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real Ts = pin->GetReal("problem", "Ts");
  Real Qs = pin->GetReal("problem", "Qs");
  Real drho = pin->GetReal("problem", "drho");
  Real zmin = pin->GetReal("problem", "zmin");
  Real zmax = pin->GetReal("problem", "zmax");

  // construct a 1D pseudo-moist adiabat with given relative humidity
  Real x1min = pmy_mesh->mesh_size.x1min;
  Real x1max = pmy_mesh->mesh_size.x1max;
  int nx1 = 2*pmy_mesh->mesh_size.nx1 + 1;
  Real dz = (x1max - x1min)/(nx1 - 1);
  Real **w1, *z1;
  NewCArray(w1, nx1, NHYDRO+2*NVAPOR);
  z1 = new Real [nx1];

  std::fill(w1[0], w1[0] + NHYDRO+2*NVAPOR, 0.);
  w1[0][iH2O] = Qs;
  Real Ps = p0;
  pthermo->ConstructAtmosphere(w1, Ts, Ps, grav, dz, nx1, Adiabat::reversible);
  z1[0] = x1min;
  for (int i = 1; i < nx1; ++i)
    z1[i] = z1[i-1] + dz;

  // setup initial condition
  Particles *pp = ppart->FindParticle("H2Op");
  for (int i = is; i <= ie; ++i) {
    Real buf[NHYDRO+2*NVAPOR];
    interpn(buf, &pcoord->x1v(i), *w1, z1, &nx1, 1, NHYDRO+2*NVAPOR);
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j) {
        for (int n = 0; n < NHYDRO; ++n) {
          phydro->w(n,k,j,i) = buf[n];
        }
        Real x1 = pcoord->x1v(i) - 8.E3;
        Real x2 = pcoord->x2v(j) - 2.E3;
        if (sqrt(x1*x1 + x2*x2) < 400.)
          pp->c(0,k,j,i) = buf[NHYDRO] + buf[NHYDRO+1];
        x1 = pcoord->x1v(i) - 5.E3;
        x2 = pcoord->x2v(j) - 5.E3;
        phydro->w(IM1,k,j,i) = -omega*x2;
        phydro->w(IM2,k,j,i) = omega*x1;
      }
    //Real cloud  = buf[NHYDRO]+buf[NHYDRO+NVAPOR];
    //std::cout << (buf[IDN]*buf[iH2O]+cloud)/(buf[IDN]+cloud) << " ";
    //std::cout << std::endl;
  }

  /* add additional precipitation 
  for (int i = is; i <= ie; ++i) {
    Real c1 = pcoord->x1f(i), c2 = pcoord->x1f(i+1);
    Real dc = c2 - c1;
    if ((c2 > zmin) && (c1 < zmax)) {
      Real f1 = (c2 - zmin)/dc;
      Real f2 = (zmax - c1)/dc;
      phydro->w(iH2Op,i) = std::min(phydro->w(iH2Oc,i), min(1.,f1,f2)*drho/phydro->w(IDN,i));
      phydro->w(iH2Oc,i) -= phydro->w(iH2Op,i);
    }
  }*/

  ppart->Initialize();
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
  FreeCArray(w1);
  delete[] z1;
}
