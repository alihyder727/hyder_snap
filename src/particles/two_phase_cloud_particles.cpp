/** @file two_phase_cloud_particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 13:38:00 PDT
 * @bug No known bugs.
 */

#include "../mesh/mesh.hpp"
#include "particles.hpp"

TwoPhaseCloudParticles::TwoPhaseCloudParticles(MeshBlock *pmb, ParameterInput *pin):
  Particles(pmb, pin, "2pcp")
{
  int nc1 = lengths_[2], nc2 = lengths_[1], nc3 = lengths_[0];
  c.NewAthenaArray(2*NVAPOR, nc3, nc2, nc1);
  dc.NewAthenaArray(2*NVAPOR, nc3, nc2, nc1);

  max_number_ = pin->GetOrAddInteger("particles", "2pcp.max_number", 1<<20);
  seeds_per_cell_ = pin->GetOrAddInteger("particles", "2pcp.seeds_per_cell", 1<<6);
}

