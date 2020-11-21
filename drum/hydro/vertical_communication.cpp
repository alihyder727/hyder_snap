// Athena++ headers
#include "hydro.hpp"
#include "vertical_communication.hpp"
#include "../utils/utils.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

#define MAX_DATA_SIZE 64

VerticalCommunication::VerticalCommunication(Hydro *phydro) :
    pmy_hydro_(phydro)
{
  MeshBlock *pmb = phydro->pmy_block;
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  buffer_ = new Real [MAX_DATA_SIZE];
  NewCArray(implicit_data_, nc3, nc2, nc1, MAX_DATA_SIZE);
#ifdef MPI_PARALLEL
  NewCArray(req_send_bot_data_, nc3, nc2);
  NewCArray(req_send_top_data_, nc3, nc2);
#endif
}

VerticalCommunication::~VerticalCommunication() {
  delete[] buffer_;
  FreeCArray(implicit_data_);
#ifdef MPI_PARALLEL
  FreeCArray(req_send_bot_data_);
  FreeCArray(req_send_top_data_);
#endif
}

void VerticalCommunication::FindNeighbors() {
  // find top and bot neighbor
  has_top_neighbor = false;
  has_bot_neighbor = false;
  for (int n = 0; n < pmy_hydro_->pmy_block->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmy_hydro_->pmy_block->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      bblock = nb;
      has_bot_neighbor = true;
    } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      tblock = nb;
      has_top_neighbor = true;
    }
  }
}

void VerticalCommunication::WaitSendTopImplicit(int kl, int ku, int jl, int ju) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      MPI_Wait(&req_send_top_data_[k][j], &status);
#endif
}

void VerticalCommunication::WaitSendBotImplicit(int kl, int ku, int jl, int ju) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      MPI_Wait(&req_send_bot_data_[k][j], &status);
#endif
}

int VerticalCommunication::CreateMPITag(int lid, int bufid, int phys) {
  return (lid<<17) | (bufid<<11) | phys;
}

#undef MAX_DATA_SIZE
