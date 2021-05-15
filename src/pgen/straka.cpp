/*======================================================================================
 * Athena++/Atmosphere Example Program
 *
 * Author: Cheng Li, University of Michigan, 2021
 * Reference: Straka et al., 1993
 *======================================================================================
 */

// @sect3{Include files}

// These are Athena++ headers files
// Athena++ is able to run both single precision (float) and double precision (double) applications
// A macro (Real) is defined to indicate the precision and it is define here.
#include "../athena.hpp"

// This file defines the basic multi-dimension array that stores fluid dynamic data.
#include "../athena_arrays.hpp"

// The input file and parameters are treated by the ParameterInput class and is define
// here.
#include "../parameter_input.hpp"

// Boundary condition is taken care of by the BoundaryValues class defined in the file
// below.
#include "../bvals/bvals.hpp"

// Coordinate related information are stored in the Coordinates class.
#include "../coordinates/coordinates.hpp"

// This file describes the equation of state.
#include "../eos/eos.hpp"

// Magnetohydrodynamics are treated here. 
// Since this is hydro-only application, this file is only a place holder.
#include "../field/field.hpp"

// Hydrodynamic fields are stored in the Hydro class
#include "../hydro/hydro.hpp"

// dynamic fields are solved on a mesh and 
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"
#include "../thermodynamics/thermodynamics.hpp"

// This is the math library. It is better to put math-related functions in a separate
// namespace
namespace math {
  #include "../math/core.h"
}

// Defines global variables
Real K, p0, cp, Rd;

// @sect3{User output variables}
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(2);
  SetUserOutputVariableName(0, "temp");
  SetUserOutputVariableName(1, "theta");
}

// @sect3{Calculate user output variables}
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
{
  for (int j = js; j <= je; ++j)
    for (int i = is; i <= ie; ++i) {
      Real prim[NHYDRO];
      for (int n = 0; n < NHYDRO; ++n)
        prim[n] = phydro->w(n,j,i);
      user_out_var(0,j,i) = pthermo->Temp(phydro->w.at(j,i));
      user_out_var(1,j,i) = pthermo->Theta(phydro->w.at(j,i), p0);
    }
}

// @sect3{Forcing function}
void Diffusion(MeshBlock *pmb, Real const time, Real const dt,
  AthenaArray<Real> const& w, AthenaArray<Real> const& bcc, AthenaArray<Real> &u)
{
  Real dx = pmb->pcoord->dx1f(0);
  Thermodynamics *pthermo = pmb->pthermo;
  for (int j = pmb->js; j <= pmb->je; ++j)
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      Real prim[NHYDRO];
      Real temp = pthermo->Temp(w.at(j,i));
      Real theta = pthermo->Theta(w.at(j,i), p0);
      Real theta_ip1_j = pthermo->Theta(w.at(j+1,i), p0);
      Real theta_im1_j = pthermo->Theta(w.at(j-1,i), p0);
      Real theta_i_jp1 = pthermo->Theta(w.at(j,i+1), p0);
      Real theta_i_jm1 = pthermo->Theta(w.at(j,i-1), p0);

      u(IM1,j,i) += dt*K*w(IDN,j,i)/(dx*dx)*(
        w(IV1,j,i-1) + w(IV1,j,i+1) + w(IV1,j-1,i) + w(IV1,j+1,i) - 4.*w(IV1,j,i));
      u(IM2,j,i) += dt*K*w(IDN,j,i)/(dx*dx)*(
        w(IV2,j,i-1) + w(IV2,j,i+1) + w(IV2,j-1,i) + w(IV2,j+1,i) - 4.*w(IV2,j,i));
      u(IEN,j,i) += dt*K*w(IDN,j,i)/(dx*dx)*cp*temp/theta*(theta_ip1_j + theta_im1_j +
        theta_i_jp1 + theta_i_jm1 - 4.*theta);
    }
}

// @sect3{Set program varialbes here}
void Mesh::InitUserMeshData(ParameterInput *pin)
{
  // forcing parameters
  Real gamma = pin->GetReal("hydro", "gamma");
  K  = pin->GetReal("problem", "K");
  p0 = pin->GetReal("problem", "p0");
  Rd = pin->GetReal("thermodynamics", "Rd");
  cp = gamma/(gamma - 1.)*Rd;

  // forcing function
  EnrollUserExplicitSourceFunction(Diffusion);
}

// @sect3{Main problem generator}
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  Real grav = -phydro->hsrc.GetG1();
  Real Ts = pin->GetReal("problem", "Ts");
  Real xc = pin->GetReal("problem", "xc");
  Real xr = pin->GetReal("problem", "xr");
  Real zc = pin->GetReal("problem", "zc");
  Real zr = pin->GetReal("problem", "zr");
  Real dT = pin->GetReal("problem", "dT");

  for (int j = js; j <= je; ++j) {
    for (int i = is; i <= ie; ++i) {
      Real x1 = pcoord->x1v(i);
      Real x2 = pcoord->x2v(j);
      Real L = sqrt(math::sqr((x2 - xc)/xr) + math::sqr((x1 - zc)/zr));
      Real temp = Ts - grav*x1/cp;
      phydro->w(IPR,j,i) = p0*pow(temp/Ts, cp/Rd);
      if (L <= 1.)
        temp += dT*(cos(M_PI*L) + 1.)/2.;
      phydro->w(IDN,j,i) = phydro->w(IPR,j,i)/(Rd*temp);
      phydro->w(IV1,j,i) = phydro->w(IV2,j,i) = 0.;
    }
  }

  if (pbval->block_bcs[inner_x1] == BoundaryFlag::outflow) {
    for (int j = js-1; j <= je+1; ++j)
      for (int i = 1; i <= NGHOST; ++i) {
        Real temp = Ts - grav*pcoord->x1f(is)/cp;
        phydro->w(IPR,j,is-i) = p0*pow(temp/Ts, cp/Rd);
        phydro->w(IDN,j,is-i) = phydro->w(IPR,j,is-i)/(Rd*temp);
        phydro->w(IVX,j,is-i) = phydro->w(IVY,j,is-i) = 0.;
      }
  }

  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
