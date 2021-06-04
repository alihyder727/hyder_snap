/** @file two_phase_cloud_particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Monday May 31, 2021 13:38:00 PDT
 * @bug No known bugs.
 */

#include "../mesh/mesh.hpp"
#include "particles.hpp"

TwoPhaseCloudParticles::TwoPhaseCloudParticles(
  MeshBlock *pmb, ParameterInput *pin, std::string name):
  Particles(pmb, pin, name)
{
  int nc1 = dims_[2], nc2 = dims_[1], nc3 = dims_[0];
  c.NewAthenaArray(2, nc3, nc2, nc1);
  dc.NewAthenaArray(2, nc3, nc2, nc1);
  cnames_.resize(2);
  cnames_[0] = "liquid";
  cnames_[1] = "solid";

  max_number_ = pin->GetOrAddInteger("particles", name + ".max_number", 1<<20);
  seeds_per_cell_ = pin->GetOrAddInteger("particles", name + ".seeds_per_cell", 1<<6);
}
