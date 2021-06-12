/** @file simple_cloud_particles.cpp
 * @brief
 *
 * @author Cheng Li (chengcli@umich.edu)
 * @date Wednesday Jun 09, 2021 18:33:03 PDT
 * @bug No known bugs.
 */

#include "../mesh/mesh.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "particles.hpp"

SimpleCloudParticles::SimpleCloudParticles(
  MeshBlock *pmb, ParameterInput *pin, std::string name):
  Particles(pmb, pin, name, 2)
{
  cnames_.resize(2);
  cnames_[0] = "cloud";
  cnames_[1] = "rain";

  seeds_per_cell_ = pin->GetOrAddInteger("particles", name + ".seeds_per_cell", 10);

  Real mu = pin->GetReal("particles", name + ".mu");
  Real cc = pin->GetReal("particles", name + ".cc");

  mu_.push_back(mu);
  mu_.push_back(mu);

  cc_.push_back(cc);
  cc_.push_back(cc);
}
