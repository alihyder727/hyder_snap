/** @file material_point.hpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 14:23:40 PDT
 * @bug No known bugs.
 */

#ifndef MATERIAL_POINT_HPP
#define MATERIAL_POINT_HPP

// C/C++ headers
#include <iosfwd>

// Athena++ headers
#include "../athena.hpp"

struct MaterialPoint {
  //! alive flag
  int alive;
  //! cell id
  int ci, cj, ck;
  //! category id
  int ct;

  #if NINT_PARTICLE_DATA > 0
    int  ii[NINT_PARTICLE_DATA];
  #endif

  //! mass in [kg]
  Real mass;
  //! positions
  Real x1, x2, x3;
  //! velocities
  Real v1, v2, v3;
  //! accelerations
  Real a1, a2, a3;

  #if NREAL_PARTICLE_DATA > 0
    Real rr[NREAL_PARTICLE_DATA];
  #endif

  MaterialPoint();
  ~MaterialPoint();
  MaterialPoint(MaterialPoint const& other);
  MaterialPoint& operator=(MaterialPoint const& other);
};

std::ostream& operator<<(std::ostream &os, MaterialPoint const& mp);

#endif /* end of include guard MATERIAL_POINT_HPP */

