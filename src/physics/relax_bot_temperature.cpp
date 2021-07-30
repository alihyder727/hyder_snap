#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "physics.hpp"

TaskStatus Physics::RelaxBotTemperature(AthenaArray<Real> &du,
  AthenaArray<Real> const& w, Real time, Real dt)
{
  MeshBlock *pmb = pmy_block;
  Thermodynamics *pthermo = pmb->pthermo;

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  Real Rd = pthermo->GetRd();
  Real gamma = pmb->peos->GetGamma();

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      Real cv = pthermo->GetMeanCv(w.at(k,j,ie));
      Real tem = pthermo->GetTemp(w.at(k,j,js));
      du(IEN,k,j,is) += dt/tau_Tbot_*(tem_bot_(k,j) - tem)*cv;
    }

  return TaskStatus::success;
}
