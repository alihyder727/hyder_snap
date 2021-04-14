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
  int nc1, nc2, nc3;
  if (dir == X1DIR) {
    nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  } else if (dir == X2DIR) {
    nc1 = pmb->ncells2, nc2 = pmb->ncells3, nc3 = pmb->ncells1;
  } else { // X3DIR
    nc1 = pmb->ncells3, nc2 = pmb->ncells1, nc3 = pmb->ncells2;
  }

  du_.NewAthenaArray(NHYDRO, nc3, nc2, nc1);
  du_.ZeroClear();
  buffer_ = new Real [MAX_DATA_SIZE];
  usend_top_ = new Real [NHYDRO*nc3*nc2];
  urecv_bot_ = new Real [NHYDRO*nc3*nc2];
  usend_bot_ = new Real [NHYDRO*nc3*nc2];
  urecv_top_ = new Real [NHYDRO*nc3*nc2];
  NewCArray(coefficients_, nc3, nc2, nc1, MAX_DATA_SIZE);

#ifdef MPI_PARALLEL
  NewCArray(req_send_bot_data_, nc3, nc2);
  NewCArray(req_send_top_data_, nc3, nc2);
#endif
}

ImplicitSolver::~ImplicitSolver() {
  delete[] buffer_;
  delete[] usend_top_;
  delete[] urecv_bot_;
  delete[] usend_bot_;
  delete[] urecv_top_;
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

void ImplicitSolver::SynchronizeConserved(AthenaArray<Real> const& du,
  int kl, int ku, int jl, int ju, int is, int ie) {
  MeshBlock *pmb = pmy_hydro->pmy_block;

  if (has_bot_neighbor) {
    int sbot = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int n = 0; n < NHYDRO; ++n)
          usend_bot_[sbot++] = du(n,k,j,is);
    if (bblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 17);
      MPI_Isend(usend_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_sync_bot_);
      tag = CreateMPITag(pmb->gid, bblock.snb.gid, 19);
      MPI_Irecv(urecv_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_recv_sync_bot_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      if (mydir == X1DIR)
        std::memcpy(pbl->phydro->pimp1->urecv_top_, usend_bot_, sbot*sizeof(Real));
      else if (mydir == X2DIR)
        std::memcpy(pbl->phydro->pimp2->urecv_top_, usend_bot_, sbot*sizeof(Real));
      else
        std::memcpy(pbl->phydro->pimp3->urecv_top_, usend_bot_, sbot*sizeof(Real));
    }
  }

  if (has_top_neighbor) {
    int stop = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int n = 0; n < NHYDRO; ++n)
          usend_top_[stop++] = du(n,k,j,ie);
    if (tblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(tblock.snb.gid, pmb->gid, 19);
      MPI_Isend(usend_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_sync_top_);
      tag = CreateMPITag(pmb->gid, tblock.snb.gid, 17);
      MPI_Irecv(urecv_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_recv_sync_top_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      if (mydir == X1DIR)
        std::memcpy(pbl->phydro->pimp1->urecv_bot_, usend_top_, stop*sizeof(Real));
      else if (mydir == X2DIR)
        std::memcpy(pbl->phydro->pimp2->urecv_bot_, usend_top_, stop*sizeof(Real));
      else
        std::memcpy(pbl->phydro->pimp3->urecv_bot_, usend_top_, stop*sizeof(Real));
    }
  }
}

void ImplicitSolver::WaitToFinishSync(int kl, int ku, int jl, int ju, int is, int ie) {
#ifdef MPI_PARALLEL
  MPI_Status status;
  if (has_bot_neighbor && (bblock.snb.rank != Globals::my_rank)) {
    MPI_Wait(&req_send_sync_bot_, &status);
    MPI_Wait(&req_recv_sync_bot_, &status);
  }
#endif

  if (has_bot_neighbor) {
    int p = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int n = 0; n < NHYDRO; ++n)
          du_(n,k,j,is-1) = urecv_bot_[p++];
  } else {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int n = 0; n < NHYDRO; ++n)
          du_(n,k,j,is-1) = du_(n,k,j,is);
  }
  
#ifdef MPI_PARALLEL
  if (has_top_neighbor && (tblock.snb.rank != Globals::my_rank)) {
    MPI_Wait(&req_send_sync_top_, &status);
    MPI_Wait(&req_recv_sync_top_, &status);
  }
#endif

  if (has_top_neighbor) {
    int p = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int n = 0; n < NHYDRO; ++n)
          du_(n,k,j,ie+1) = urecv_top_[p++];
  } else {
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int n = 0; n < NHYDRO; ++n)
          du_(n,k,j,ie+1) = du_(n,k,j,ie);
  }
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

void ImplicitSolver::WaitToFinishSend(int kl, int ku, int jl, int ju)
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
