// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../thermodynamics/thermodynamics.hpp"
#include "../../../utils/utils.hpp"
#include "../../../hydro/hydro.hpp"
#include "bvals_hydro.hpp"

void HydroBoundaryVariable::OutflowInnerX1(
    Real time, Real dt, int il, int jl, int ju, int kl, int ku, int ngh)
{
  // extend velocities
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_cc)(IVX,k,j,il-i) = (*var_cc)(IVX,k,j,il)/pow(2,i-1);
        (*var_cc)(IVY,k,j,il-i) = (*var_cc)(IVY,k,j,il)/pow(2,i-1);
        (*var_cc)(IVZ,k,j,il-i) = (*var_cc)(IVZ,k,j,il)/pow(2,i-1);
      }

  // conforms to the case without gravity
  if (pmy_block_->phydro->hsrc.GetG1() == 0.) {
    for (int n = 0; n <= NVAPOR; ++n) {
      for (int k=kl; k<=ku; ++k)
        for (int j=jl; j<=ju; ++j)
#pragma omp simd
          for (int i=1; i<=ngh; ++i)
            (*var_cc)(n,k,j,il-i) = (*var_cc)(n,k,j,il);
    }
    for (int k=kl; k<=ku; ++k)
      for (int j=jl; j<=ju; ++j)
#pragma omp simd
        for (int i=1; i<=ngh; ++i)
          (*var_cc)(IPR,k,j,il-i) = (*var_cc)(IPR,k,j,il);
  }
}

void HydroBoundaryVariable::OutflowOuterX1(
    Real time, Real dt, int iu, int jl, int ju, int kl, int ku, int ngh)
{
  // extend velocities
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        (*var_cc)(IVX,k,j,iu+i) = (*var_cc)(IVX,k,j,iu)/pow(2,i-1);
        (*var_cc)(IVY,k,j,iu+i) = (*var_cc)(IVY,k,j,iu)/pow(2,i-1);
        (*var_cc)(IVZ,k,j,iu+i) = (*var_cc)(IVZ,k,j,iu)/pow(2,i-1);
      }

  // conforms to the case without gravity
  if (pmy_block_->phydro->hsrc.GetG1() == 0.) {
    for (int n = 0; n <= NVAPOR; ++n) {
      for (int k=kl; k<=ku; ++k)
        for (int j=jl; j<=ju; ++j)
#pragma omp simd
          for (int i=1; i<=ngh; ++i)
            (*var_cc)(n,k,j,iu+i) = (*var_cc)(n,k,j,iu);
    }
    for (int k=kl; k<=ku; ++k)
      for (int j=jl; j<=ju; ++j)
#pragma omp simd
        for (int i=1; i<=ngh; ++i)
          (*var_cc)(IPR,k,j,iu+i) = (*var_cc)(IPR,k,j,iu);
  }
}
