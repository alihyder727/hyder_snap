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
  Particles(pmb, pin, name, 2)
{
  cnames_.resize(2);
  cnames_[0] = "liquid";
  cnames_[1] = "solid";

  max_number_ = pin->GetOrAddInteger("particles", name + ".max_number", 1<<20);
  seeds_per_cell_ = pin->GetOrAddInteger("particles", name + ".seeds_per_cell", 1<<6);
}

void TwoPhaseCloudParticles::ExchangeHydro(std::vector<MaterialPoint> &mp,
  AthenaArray<Real> &du, AthenaArray<Real> const &w)
{
  Particles::ExchangeHydro(mp, du, w);
  for (std::vector<MaterialPoint>::iterator it = mp.begin(); it != mp.end(); ++it) {
    it->v1 += -10.;
  }
}
