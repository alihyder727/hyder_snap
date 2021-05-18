#include "diagnostics.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../coordinates/coordinates.hpp"
#include "../reconstruct/interpolation.hpp"

Buoyancy::Buoyancy(MeshBlock *pmb) : Diagnostics(pmb, "b")
{
  type = "SCALARS";
  data.NewAthenaArray(ncells3_,ncells2_,ncells1_);
  wl_.NewAthenaArray(NHYDRO,ncells3_,ncells2_,ncells1_+1);
  wr_.NewAthenaArray(NHYDRO,ncells3_,ncells2_,ncells1_+1);
  grav_ = pmb->phydro->hsrc.GetG1();
}

Buoyancy::~Buoyancy()
{
  data.DeleteAthenaArray();
  wl_.DeleteAthenaArray();
  wr_.DeleteAthenaArray();
}

void Buoyancy::Finalize(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block_;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie+1; ++i) {
        wl_(IPR,k,j,i+1) = interp_weno5(w(IPR,k,j,i+2),w(IPR,k,j,i+1),
                                        w(IPR,k,j,i),w(IPR,k,j,i-1),w(IPR,k,j,i-2));
        wr_(IPR,k,j,i) = interp_weno5(w(IPR,k,j,i-2),w(IPR,k,j,i-1),
                                      w(IPR,k,j,i),w(IPR,k,j,i+1),w(IPR,k,j,i+2));
      }

  // buoyancy
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real dx = pmb->pcoord->dx1f(i);
        data(k,j,i) = -(wl_(IPR,k,j,i+1) - wr_(IPR,k,j,i))/(dx*w(IDN,k,j,i)) + grav_;
      }
}
