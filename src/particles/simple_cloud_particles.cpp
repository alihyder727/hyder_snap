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
  Thermodynamics *pthermo = pmb->pthermo;
  cnames_.resize(2);
  cnames_[0] = "cloud";
  cnames_[1] = "rain";

  seeds_per_cell_ = pin->GetOrAddInteger("particles", name + ".seeds_per_cell", 1<<6);

  int ic = pin->GetInteger("particles", name + ".index");

  mu_.push_back(pthermo->GetMassRatio(ic)*Thermodynamics::Rgas/pthermo->GetRd());
  mu_.push_back(pthermo->GetMassRatio(ic)*Thermodynamics::Rgas/pthermo->GetRd());

  cc_.push_back(pthermo->GetCp(ic));
  cc_.push_back(pthermo->GetCp(ic));
}
