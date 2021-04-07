#include "../mesh/mesh.hpp"
#include "../eos/eos.hpp"
#include "../hydro/hydro.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"
#include "physics.hpp"

TaskStatus Physics::RelaxBotTemperature(AthenaArray<Real> &u,
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
      Real v1 = u(IM1,k,j,is);
      Real v2 = u(IM2,k,j,is);
      Real v3 = u(IM3,k,j,is);
      Real pres = u(IDN,k,j,is)*Rd*tem_bot_(k,j);
      Real KE = 0.5*(v1*v1+v2*v2+v3*v3)/u(IDN,k,j,is);
      Real LE = 0., fsig = 1., feps = 1.;

      // clouds
      for (int n = 1 + NVAPOR; n < NMASS; ++n)
        LE += -pthermo->GetLatent(n)*u(n,k,j,is);

      // vapors
      for (int n = 1; n <= NVAPOR; ++n) {
        LE += -pthermo->GetLatent(n)*u(n,k,j,is);
        fsig += u(n,k,j,is)/u(IDN,k,j,is)*pthermo->GetCvRatio(n);
        feps += u(n,k,j,is)/u(IDN,k,j,is)/pthermo->GetMassRatio(n);
        pres += u(n,k,j,is)*Rd/pthermo->GetMassRatio(n)*tem_bot_(k,j);
      }

      Real u0 = pres/(gamma - 1.)*fsig/feps + KE + LE;
      u(IEN,k,j,is) = (u(IEN,k,j,is) + dt/tau_Tbot_*u0)/(1. + dt/tau_Tbot_);
    }

  return TaskStatus::success;
}
