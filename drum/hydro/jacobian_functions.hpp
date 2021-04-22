#ifndef JACOBIAN_FUNCTIONS_HPP
#define JACOBIAN_FUNCTIONS_HPP

// Athena++ header
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"
#include "hydro.hpp"

template<typename T>
void JacobianGravityCoriolis(T &jac, int k, int j, int i,
  Hydro *phydro, CoordinateDirection dir) {
  Real omega1 = 0, omega2 = 0, omega3 = 0, theta, phi;
  Real omegax = phydro->hsrc.GetOmegaX();
  Real omegay = phydro->hsrc.GetOmegaY();
  Real omegaz = phydro->hsrc.GetOmegaZ();
  Real grav = phydro->hsrc.GetG1();

  MeshBlock *pmb = phydro->pmy_block;

  if (COORDINATE_SYSTEM == "cartesian") {
    omega1 = omegax;
    omega2 = omegay;
    omega3 = omegaz;
  } else if (COORDINATE_SYSTEM == "cylindrical") {
    theta = pmb->pcoord->x2v(j);

    omega1 = cos(theta)*omegax + sin(theta)*omegay;
    omega2 = -sin(theta)*omegax + cos(theta)*omegay;
    omega3 = omegaz;
  } else if (COORDINATE_SYSTEM == "spherical_polar") {
    theta = pmb->pcoord->x2v(j);
    phi = pmb->pcoord->x3v(k);

    omega1 = sin(theta)*cos(phi)*omegax + sin(theta)*sin(phi)*omegay + cos(theta)*omegaz;
    omega2 = cos(theta)*cos(phi)*omegax + cos(theta)*sin(phi)*omegay - sin(theta)*omegaz;
    omega3 = -sin(phi)*omegax + cos(phi)*omegay;
  }

  jac.setZero();
  int idn = 0, ivx = 1, ivy = 2, ivz = 3, ien = 4;
  if (dir == X1DIR) {
    jac(ivx,idn) = grav;
    jac(ien,ivx) = grav;
  } else if (dir == X2DIR) {
    Real tmp = omega1;
    omega1 = omega2;
    omega2 = omega3;
    omega3 = tmp;
  } else { // X3DIR
    Real tmp = omega3;
    omega3 = omega2;
    omega2 = omega1;
    omega1 = tmp;
  }

  jac(ivx,ivy) = 2.*omega3;
  jac(ivx,ivz) = -2*omega2;
  jac(ivy,ivx) = -2*omega3;
  jac(ivy,ivz) = 2*omega1;
  jac(ivz,ivx) = 2*omega2;
  jac(ivz,ivy) = -2*omega1;
}

#endif
