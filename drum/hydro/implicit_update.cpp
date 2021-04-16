//! \file implicit_update.cpp
//  \brief final step of full implicit solver

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "hydro.hpp"
#include "implicit/implicit_solver.hpp"

void Hydro::ImplicitUpdate(AthenaArray<Real> &du)
{
  MeshBlock *pmb = pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  if (implicit_flag == 1) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,k,j,i) = pimp1->du_(IDN,k,j,i);
          du(IVX,k,j,i) = pimp1->du_(IVX,k,j,i);
          du(IEN,k,j,i) = pimp1->du_(IEN,k,j,i);
        }
  } else if (implicit_flag == 2) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,k,j,i) = pimp2->du_(IDN,i,k,j);
          du(IVY,k,j,i) = pimp2->du_(IVY,i,k,j);
          du(IEN,k,j,i) = pimp2->du_(IEN,i,k,j);
        }
  } else if (implicit_flag == 3) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,k,j,i) = pimp3->du_(IDN,j,i,k);
          du(IVZ,k,j,i) = pimp3->du_(IVZ,j,i,k);
          du(IEN,k,j,i) = pimp3->du_(IEN,j,i,k);
        }
  } else if (implicit_flag == 4) {
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          du(IDN,k,j,i) = pimp1->du_(IDN,k,j,i) + pimp2->du_(IDN,i,k,j)
            + pimp3->du_(IDN,j,i,k) - 2*du(IDN,k,j,i);
          du(IVX,k,j,i) = pimp1->du_(IVX,k,j,i);
          du(IVY,k,j,i) = pimp2->du_(IVY,i,k,j);
          du(IVZ,k,j,i) = pimp3->du_(IVZ,j,i,k);
          du(IEN,k,j,i) = pimp1->du_(IEN,k,j,i) + pimp2->du_(IEN,i,k,j)
            + pimp3->du_(IEN,j,i,k) - 2*du(IEN,k,j,i);
        }
  }
}
