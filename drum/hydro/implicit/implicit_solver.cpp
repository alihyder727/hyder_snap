// C/C++ headers
#include <string>
#include <functional>

// Athena++ headers
#include "../hydro.hpp"
#include "../../utils/utils.hpp"
#include "implicit_solver.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

#define MAX_DATA_SIZE 64

ImplicitSolver::ImplicitSolver(Hydro *phydro, CoordinateDirection dir):
    pmy_hydro(phydro), mydir(dir), has_bot_neighbor(false), has_top_neighbor(false)
{
  MeshBlock *pmb = phydro->pmy_block;
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  buffer_ = new Real [MAX_DATA_SIZE];
  if (dir == X1DIR) {
    NewCArray(coefficients_, nc3, nc2, nc1, MAX_DATA_SIZE);
#ifdef MPI_PARALLEL
    NewCArray(req_send_bot_data_, nc3, nc2);
    NewCArray(req_send_top_data_, nc3, nc2);
    du_.NewAthenaArray(NHYDRO, nc3, nc2, nc1);
    du_.ZeroClear();
#endif
  } else if (dir == X2DIR) {
    NewCArray(coefficients_, nc1, nc3, nc2, MAX_DATA_SIZE);
#ifdef MPI_PARALLEL
    NewCArray(req_send_bot_data_, nc1, nc3);
    NewCArray(req_send_top_data_, nc1, nc3);
    du_.NewAthenaArray(NHYDRO, nc1, nc3, nc2);
    du_.ZeroClear();
#endif
  } else { // X3DIR
    NewCArray(coefficients_, nc2, nc1, nc3, MAX_DATA_SIZE);
#ifdef MPI_PARALLEL
    NewCArray(req_send_bot_data_, nc2, nc1);
    NewCArray(req_send_top_data_, nc2, nc1);
    du_.NewAthenaArray(NHYDRO, nc2, nc1, nc3);
    du_.ZeroClear();
#endif
  }
}

ImplicitSolver::~ImplicitSolver() {
  delete[] buffer_;
  FreeCArray(coefficients_);
#ifdef MPI_PARALLEL
  FreeCArray(req_send_bot_data_);
  FreeCArray(req_send_top_data_);
#endif
}

void ImplicitSolver::FindNeighbors() {
  // find top and bot neighbor
  has_top_neighbor = false;
  has_bot_neighbor = false;
  for (int n = 0; n < pmy_hydro->pmy_block->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmy_hydro->pmy_block->pbval->neighbor[n];
    if (mydir == X1DIR) {
      if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
        bblock = nb;
        has_bot_neighbor = true;
      } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
        tblock = nb;
        has_top_neighbor = true;
      }
    } else if (mydir == X2DIR) {
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == -1) && (nb.ni.ox3 == 0)) {
        bblock = nb;
        has_bot_neighbor = true;
      } if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 1) && (nb.ni.ox3 == 0)) {
        tblock = nb;
        has_top_neighbor = true;
      }
    } else { // X3DIR
      if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == -1)) {
        bblock = nb;
        has_bot_neighbor = true;
      } if ((nb.ni.ox1 == 0) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 1)) {
        tblock = nb;
        has_top_neighbor = true;
      }
    }
  }
}

void ImplicitSolver::SynchronizeConserved(AthenaArray<Real> &du_) {
}

void ImplicitSolver::WaitSendTop(int kl, int ku, int jl, int ju) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      MPI_Wait(&req_send_top_data_[k][j], &status);
#endif
}

void ImplicitSolver::WaitSendBot(int kl, int ku, int jl, int ju) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      MPI_Wait(&req_send_bot_data_[k][j], &status);
#endif
}

void ImplicitSolver::WaitToFinish(int kl, int ku, int jl, int ju)
{
  if (has_top_neighbor)
    WaitSendTop(kl, ku, jl, ju);

  if (has_bot_neighbor)
    WaitSendBot(kl, ku, jl, ju);
}

int ImplicitSolver::CreateMPITag(int recvid, int sendid, int phys) {
  //return (lid<<17) | (bufid<<11) | phys;
  std::string str = std::to_string(recvid);
  str += std::to_string(sendid);
  str += std::to_string(phys);
  return std::hash<std::string>{}(str)%(1<<24);
}

#undef MAX_DATA_SIZE
