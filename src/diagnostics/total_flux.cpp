#include "diagnostics.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../thermodynamics/thermodynamics.hpp"

TotalFlux::TotalFlux(MeshBlock *pmb):Diagnostics(pmb, "flux")
{
  type = "VECTORS";
  grid = "--C";
  data.NewAthenaArray(NHYDRO,1,1,ncells1_);
}

TotalFlux::~TotalFlux()
{
  data.DeleteAthenaArray();
}

void TotalFlux::Progress(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  if (ncycle == 0)
    std::fill(data.data(), data.data() + data.GetSize(), 0.);

  // sum over horizontal grids weighted by volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k,j,is,ie,vol_);
      for (int i = is; i <= ie; ++i) {
        data(0,i) += vol_(i)*w(IDN,k,j,i)*w(IVX,k,j,i);
        for (int n = 1; n < IPR; ++n)
          data(n,i) += vol_(i)*w(n,k,j,i)*w(IVX,k,j,i);
        data(IPR,i) += vol_(i)*pthermo->GetTemp(w.at(k,j,i))*w(IVX,k,j,i);
      }
    }

  ncycle++;
}

void TotalFlux::Finalize(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block_;

  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *total_vol = new Real[ncells1_];
  std::fill(total_vol, total_vol + ncells1_, 0.);

  // calculate total volume
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pmb->pcoord->CellVolume(k,j,is,ie,vol_);
      for (int i = is; i <= ie; ++i)
        total_vol[i] += vol_(i);
    }

  // take time and spatial average
  if (ncycle > 0) {
  // sum over all ranks
#ifdef MPI_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, data.data(), data.GetSize(), MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, total_vol, ncells1_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif
    for (int n = 0; n < NHYDRO; ++n)
      for (int i = is; i <= ie; ++i)
        data(n,i) /= ncycle*total_vol[i];
  }

  // clear cycle;
  delete [] total_vol;
  ncycle = 0;
}
