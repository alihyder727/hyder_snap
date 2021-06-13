/** @file aggregate_density.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Jun 13, 2021 12:48:16 PDT
 * @bug No known bugs.
 */

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../math/interpolation.h" // locate
#include "particles.hpp"

void Particles::AggregateDensity(AthenaArray<Real> &c1, std::vector<MaterialPoint> const& mp)
{
  MeshBlock *pmb = pmy_block;
  Real loc[3];

  c1.ZeroClear();
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
    c1(q->type, ck, cj, ci) += q->rho;

    if (pcell_(q->type, ck, cj, ci) == nullptr) {
      pcell_(q->type,ck,cj,ci) = const_cast<MaterialPoint*>(&(*q));
      pcell_(q->type,ck,cj,ci)->next = nullptr;
    } else {  // insert sort from low density to high density
      MaterialPoint *pc = pcell_(q->type,ck,cj,ci);
      MaterialPoint *tmp;
      if (pc->rho > q->rho) {
        tmp = pcell_(q->type,ck,cj,ci);
        pcell_(q->type,ck,cj,ci) = const_cast<MaterialPoint*>(&(*q));
        pcell_(q->type,ck,cj,ci)->next = tmp;
      } else {
        while ((pc->next != nullptr) && (q->rho > pc->next->rho))
          pc = pc->next;
        tmp = pc->next;
        pc->next = const_cast<MaterialPoint*>(&(*q));
        pc->next->next = tmp;
      }
    }
  }
}

