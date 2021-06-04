/** @file translate_euler.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 13:38:00 PDT
 * @bug No known bugs.
 */

// C/C++ headers
#include <algorithm>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../bvals/bvals.hpp"
#include "particles.hpp"
#include "particle_buffer.hpp"
#include "material_point.hpp"

void Particles::TranslateEuler(std::vector<MaterialPoint> &mp, Real dt)
{
  MeshBlock *pmb = pmy_block;
  Mesh *pm = pmb->pmy_mesh;

  Real x1min = pmb->block_size.x1min;
  Real x1max = pmb->block_size.x1max;
  Real x2min = pmb->block_size.x2min;
  Real x2max = pmb->block_size.x2max;
  Real x3min = pmb->block_size.x3min;
  Real x3max = pmb->block_size.x3max;

  int ox1 = 0, ox2 = 0, ox3 = 0, fi1 = 0, fi2 = 0;
  std::vector<MaterialPoint>::iterator qi = mp.begin();
  std::vector<MaterialPoint>::iterator qj = mp.end();

  while (qi < qj) {
    qi->x1 += qi->v1*dt;
    qi->x2 += qi->v2*dt;
    qi->x3 += qi->v3*dt;

    // take care of reflective boundary condition
    if (pmb->pbval->block_bcs[inner_x1] == BoundaryFlag::reflect && qi->x1 < x1min) {
      qi->x1 = 2*x1min - qi->x1;
      qi->v1 = - qi->v1;
    }
    if (pmb->pbval->block_bcs[outer_x1] == BoundaryFlag::reflect && qi->x1 > x1max) {
      qi->x1 = 2*x1max - qi->x1;
      qi->v1 = - qi->v1;
    }
    ox1 = qi->x1 < x1min ? -1 : (qi->x1 > x1max ? 1 : 0);

    if (pm->f2 > 1) {
      if (pmb->pbval->block_bcs[inner_x2] == BoundaryFlag::reflect && qi->x2 < x2min) {
        qi->x2 = 2*x2min - qi->x2;
        qi->v2 = - qi->v2;
      }
      if (pmb->pbval->block_bcs[outer_x2] == BoundaryFlag::reflect && qi->x2 > x2max) {
        qi->x2 = 2*x2max - qi->x2;
        qi->v2 = - qi->v2;
      }
      ox2 = qi->x2 < x2min ? -1 : (qi->x2 > x2max ? 1 : 0);
    }

    if (pm->f3 > 1) {
      if (pmb->pbval->block_bcs[inner_x3] == BoundaryFlag::reflect && qi->x3 < x3min) {
        qi->x3 = 2*x3min - qi->x3;
        qi->v3 = - qi->v3;
      }
      if (pmb->pbval->block_bcs[outer_x3] == BoundaryFlag::reflect && qi->x3 > x3max) {
        qi->x3 = 2*x3max - qi->x3;
        qi->v3 = - qi->v3;
      }
      ox3 = qi->x3 < x3min ? -1 : (qi->x3 > x3max ? 1 : 0);
    }

    if (pm->multilevel) {
      // reserved implementation for multilevel, fi1, fi2
    }

    int bid = BoundaryBase::FindBufferID(ox1, ox2, ox3, fi1, fi2);

    if (qi->id > 0 && bid == -1) { // particle is alive and inside domain
      qi++;
    } else {  // particle deseased or moved out of the domain
      std::swap(*qi, *(qj-1));
      ppb->bufid.push_back(qi->id > 0 ? bid : -1); // Note that bufid is reversed
      qj--;
    }
  }
}
