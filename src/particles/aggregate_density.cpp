/** @file aggregate_density.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Jun 13, 2021 12:48:16 PDT
 * @bug No known bugs.
 */

// C/C++ headers
#include <iostream>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../math/interpolation.h" // locate
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "particles.hpp"

void Particles::AggregateDensity(AthenaArray<Real> &u, std::vector<MaterialPoint> &mp)
{
  MeshBlock *pmb = pmy_block;

  u.ZeroClear();
  std::fill(pcell_.data(), pcell_.data() + pcell_.GetSize(), nullptr);
  int i, j, k;

  for (std::vector<MaterialPoint>::iterator q = mp.begin(); q != mp.end(); ++q) {
    k = locate(xface_.data(), q->x3, dims_[0]+1);
    j = locate(xface_.data()+dims_[0]+1, q->x2, dims_[1]+1);
    i = locate(xface_.data()+dims_[0]+dims_[1]+2, q->x1, dims_[2]+1);

    assert(k >= pmb->ks && k <= pmb->ke);
    assert(j >= pmb->js && j <= pmb->je);
    assert(i >= pmb->is && i <= pmb->ie);

    u(q->type,k,j,i) += q->rho;

    if (pcell_(q->type,k,j,i) == nullptr) {
      pcell_(q->type,k,j,i) = &(*q);
      pcell_(q->type,k,j,i)->next = nullptr;
    } else {  // insert sort from low to high density
      MaterialPoint *tmp;
      if (pcell_(q->type,k,j,i)->rho < q->rho) {
        MaterialPoint *pc = pcell_(q->type,k,j,i);
        while (pc->next != nullptr && pc->next->rho < q->rho)
          pc = pc->next;
        tmp = pc->next;
        pc->next = &(*q);
        pc->next->next = tmp;
      } else {
        tmp = pcell_(q->type,k,j,i);
        pcell_(q->type,k,j,i) = &(*q);
        pcell_(q->type,k,j,i)->next = tmp;
      }
    }
  }

  // make a copy
  u1_ = u;
}
