#include "diagnostics.hpp"
#include "../globals.hpp"
#include "../thermodynamics/thermodynamics.hpp"
#include "../coordinates/coordinates.hpp"

TemperatureAnomaly::TemperatureAnomaly(MeshBlock *pmb):Diagnostics(pmb, "tempa")
{
  type = "SCALARS";
  data.NewAthenaArray(ncells3_,ncells2_,ncells1_);
}

TemperatureAnomaly::~TemperatureAnomaly()
{
  data.DeleteAthenaArray();
}

void TemperatureAnomaly::Finalize(AthenaArray<Real> const& w)
{
  MeshBlock *pmb = pmy_block_;
  Coordinates *pcoord = pmb->pcoord;
  Thermodynamics *pthermo = pmb->pthermo;
  int is = pmb->is, js = pmb->js, ks = pmb->ks;
  int ie = pmb->ie, je = pmb->je, ke = pmb->ke;

  Real *data_sum = new Real [2*ncells1_];
  std::fill(data_sum, data_sum + 2*ncells1_, 0.);

  Real *temp_mean = data_sum;
  Real *total_vol = data_sum + ncells1_;

  // calculate horizontal mean
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j) {
      pcoord->CellVolume(k,j,is,ie,vol_);
      for (int i = is; i <= ie; ++i) {
        temp_mean[i] += vol_(i)*pthermo->GetTemp(w.at(k,j,i));
        total_vol[i] += vol_(i);
      }
    }

#ifdef MPI_PARALLEL
  // sum over all ranks
  MPI_Allreduce(MPI_IN_PLACE, data_sum, 2*ncells1_, MPI_ATHENA_REAL, MPI_SUM, MPI_COMM_WORLD);
#endif

  // temperature anomaly
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i)
        data(k,j,i) = pthermo->GetTemp(w.at(k,j,i)) - temp_mean[i]/total_vol[i];

  delete [] data_sum;
}
