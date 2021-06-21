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
#include "../globals.hpp"

void Particles::Particulate(std::vector<MaterialPoint> &mp, AthenaArray<Real> const& c)
{
  MeshBlock *pmb = pmy_block;
  Coordinates *pco = pmb->pcoord;
  // use mp1 for particle buffer
  mp1.clear();
  for (int t = 0; t < c.GetDim4(); ++t)
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real delta_u = c(t,k,j,i) - c1_(t,k,j,i);
          // 1. count particles
          int nparts = 0;
          MaterialPoint *pc = pcell_(t,k,j,i);
          while (pc != nullptr) {
            pc = pc->next;
            nparts++;
          }

          // 2. merge particles
          int n0 = nparts;
          Real sum = 0.;
          pc = pcell_(t,k,j,i);
          for (int n = 0; n < n0 - nmax_per_cell_; ++n) {
            sum += pc->rho;
            available_ids_.push_back(pc->id);
            pc->id = -1;
            pc->rho = 0;
            pc = pc->next;
            nparts--;
          }
          pcell_(t,k,j,i) = pc;
          Real avg = sum/nmax_per_cell_;
          while (pc != nullptr) {
            pc->rho += avg;
            pc = pc->next;
          }

          //if (Globals::my_rank == 0)
          //  std::cout << "(" << k << "," << j << "," << i << ") = " << nparts << std::endl;
          // 3. add new particles to mp1
          if (delta_u > seeds_per_cell_*density_floor_) {
            avg = delta_u/seeds_per_cell_;
            int num = std::min(nmax_per_cell_ - nparts, seeds_per_cell_);
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
              mp1.push_back(p);
            }

            pc = pcell_(t,k,j,i);
            for (int n = 0; n < seeds_per_cell_ - num; ++n) {
              pc->rho += avg;
              pc = pc->next;
            }
          } else if (delta_u < 0) { // 4. remove particles
            Real avg = std::abs(delta_u)/nparts;
            pc = pcell_(t,k,j,i);
            //std::cout << "c =  " << c(t,k,j,i) << " c1 = " << c1_(t,k,j,i) << std::endl;
            //std::cout << "[" << pco->x1f(i) << "," << pco->x1f(i+1) << "]" << std::endl;
            //std::cout << "delta_u = " << delta_u << " nparts = " << nparts << std::endl;
            for (int n = 0; n < nparts; ++n) {
              //std::cout << pc->x1 << " " << pc->rho << " ";
              if (pc->rho > avg) {
                delta_u += avg;
                pc->rho -= avg;
              } else {
                delta_u += pc->rho;
                avg = std::abs(delta_u)/(nparts - n - 1);
                available_ids_.push_back(pc->id);
                pc->id = -1;
                pc->rho = 0;
              }
              //std::cout << pc->rho << std::endl;
              //std::cout << delta_u << " " << avg << std::endl;
              pc = pc->next;
            }
            std::stringstream msg;
            if (std::abs(delta_u) > density_floor_) {
              msg << "### FATAL ERROR in Particles::Particulate:" << std::endl
                  << "c = " << c(t,k,j,i) << " c1 = " << c1_(t,k,j,i) << std::endl
                  << "delta_u = " << delta_u << std::endl
                  << "density_floor_ = " << density_floor_ << std::endl;
              ATHENA_ERROR(msg);
            }
          }
        }

  // transfer from mp1 to mp;
  mp.reserve(mp.size() + mp1.size());
  mp.insert(mp.end(), mp1.begin(), mp1.end());
}

