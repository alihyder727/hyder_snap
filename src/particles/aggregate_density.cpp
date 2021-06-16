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
#include "particles.hpp"

void Particles::AggregateDensity(AthenaArray<Real> &c, std::vector<MaterialPoint> const& mp)
{
  MeshBlock *pmb = pmy_block;
  Real loc[3];

  c.ZeroClear();
  std::fill(pcell_.data(), pcell_.data() + pcell_.GetSize(), nullptr);
  int ci, cj, ck;

  for (std::vector<MaterialPoint>::const_iterator q = mp.begin(); q != mp.end(); ++q) {
    loc[0] = q->x3;
    loc[1] = q->x2;
    loc[2] = q->x1;

    if (dims_[0] > 1)
      ck = locate(coordinates_.data(), loc[0], dims_[0]);
    else ck = pmb->ks;

    if (dims_[1] > 1)
      cj = locate(coordinates_.data() + dims_[0], loc[1], dims_[1]);
    else cj = pmb->js;

    ci = locate(coordinates_.data() + dims_[0] + dims_[1], 
      loc[2], dims_[2]);

    assert(ck >= pmb->ks && ck <= pmb->ke);
    assert(cj >= pmb->js && cj <= pmb->je);
    assert(ci >= pmb->is && ci <= pmb->ie);

    c(q->type, ck, cj, ci) += q->rho;

    if (pcell_(q->type,ck,cj,ci) == nullptr) {
      pcell_(q->type,ck,cj,ci) = const_cast<MaterialPoint*>(&(*q));
      pcell_(q->type,ck,cj,ci)->next = nullptr;
    } else {  // insert sort from low to high density
      MaterialPoint *tmp;
      if (pcell_(q->type,ck,cj,ci)->rho < q->rho) {
        MaterialPoint *pc = pcell_(q->type,ck,cj,ci);
        while (pc->next != nullptr && pc->next->rho < q->rho)
          pc = pc->next;
        tmp = pc->next;
        pc->next = const_cast<MaterialPoint*>(&(*q));
        pc->next->next = tmp;
      } else {
        tmp = pcell_(q->type,ck,cj,ci);
        pcell_(q->type,ck,cj,ci) = const_cast<MaterialPoint*>(&(*q));
        pcell_(q->type,ck,cj,ci)->next = tmp;
      }
    }
  }

  // make a copy
  c1_ = c;
}

