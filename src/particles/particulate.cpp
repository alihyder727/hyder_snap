/** @file particulate.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Sunday Jun 13, 2021 12:16:32 PDT
 * @bug No known bugs.
 */

// C++ headers
#include <cassert>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../coordinates/coordinates.hpp"
#include "particles.hpp"

void Particles::Particulate(std::vector<MaterialPoint> &mp, AthenaArray<Real> const& c)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pco = pmb->pcoord;
  for (int t = 0; t < c.GetDim4(); ++t)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real delta_c = c(t,k,j,i) - c1_(t,k,j,i);
          int nparts = CountParticlesInCell(t,k,j,i);
          if (delta_c > density_floor_) {
            Real avg = delta_c/seeds_per_cell_;
            int num = std::min(std::max(0, nmax_per_cell_ - nparts), seeds_per_cell_);
            for (int n = 0; n < num; ++n) {
              MaterialPoint p;
              p.id = GetNextId();
              p.type = t;
              p.time = pmb->pmy_mesh->time;
              p.x1 = pco->x1f(i) + (1.*rand()/RAND_MAX)*pco->dx1f(i);
              p.x2 = pco->x2f(j) + (1.*rand()/RAND_MAX)*pco->dx2f(j);
              p.x3 = pco->x3f(k) + (1.*rand()/RAND_MAX)*pco->dx3f(k);
              p.v1 = 0.;
              p.v2 = 0.;
              p.v3 = 0.;
              p.rho = avg;
              mp.push_back(p);
            }

            MaterialPoint *pc = pcell_(t,k,j,i);
            for (int n = 0; n < seeds_per_cell_ - num; ++n) {
              pc->rho += avg;
              pc = pc->next;
            }
          } else if (delta_c < -density_floor_) {
            Real avg = std::abs(delta_c)/nparts;
            MaterialPoint *pc = pcell_(t,k,j,i);
            //std::cout << "c =  " << c(t,k,j,i) << " c1 = " << c1_(t,k,j,i) << std::endl;
            //std::cout << "[" << pco->x1f(i) << "," << pco->x1f(i+1) << "]" << std::endl;
            //std::cout << "delta_c = " << delta_c << " nparts = " << nparts << std::endl;
            while (pc != nullptr) {
              //std::cout << pc->x1 << " " << pc->rho << " ";
              nparts--;
              if (pc->rho > avg) {
                delta_c += avg;
                pc->rho -= avg;
              } else {
                delta_c += pc->rho;
                avg = std::abs(delta_c)/nparts;
                available_ids_.push_back(pc->id);
                pc->id = -1;
                pc->rho = 0;
              }
              //std::cout << pc->rho << std::endl;
              //std::cout << delta_c << " " << avg << std::endl;
              pc = pc->next;
            }
            assert(nparts == 0);
            std::stringstream msg;
            if (std::abs(delta_c) > density_floor_) {
              msg << "### FATAL ERROR in Particles::Particulate:" << std::endl
                  << "delta_c = " << delta_c << std::endl
                  << "density_floor_ = " << density_floor_ << std::endl;
              ATHENA_ERROR(msg);
            }
          }
        }
}

