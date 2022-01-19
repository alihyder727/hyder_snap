// C/C++ headers
#include <string>
#include <functional>

// Athena++ headers
#include "../hydro.hpp"
#include "../../debugger/debugger.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/utils.hpp"
#include "../../globals.hpp"
#include "implicit_solver.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

#define MAX_DATA_SIZE 25

ImplicitSolver::ImplicitSolver(Hydro *phydro, int n3max, int n2max):
    pmy_hydro(phydro), has_bot_neighbor(false), has_top_neighbor(false),
    first_block(true), last_block(true), periodic_boundary(false),
    pole_at_bot(false), pole_at_top(false)
{
  phydro->pmy_block->pdebug->Enter("ImplicitSolver");
  MeshBlock *pmb = phydro->pmy_block;
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  du_.NewAthenaArray(NHYDRO, nc3, nc2, nc1);
  coefficients_.NewAthenaArray(nc3, nc2, nc1, 4*MAX_DATA_SIZE);
  NewCArray(buffer_, n3max, n2max, 7*MAX_DATA_SIZE);

#ifdef MPI_PARALLEL
  NewCArray(req_send_data1_, n3max, n2max);
  NewCArray(req_send_data2_, n3max, n2max);
  NewCArray(req_send_data6_, n3max, n2max);
  NewCArray(req_send_data7_, n3max, n2max);
#endif

  int idn = 0, ivx = 1, ivy = 2, ivz = 3, ien = 4;
  p2_.setZero();
  p2_(idn,idn) = 1.;
  p2_(ivx,ivy) = 1.;
  p2_(ivy,ivz) = 1.;
  p2_(ivz,ivx) = 1.;
  p2_(ien,ien) = 1.;

  p3_.setZero();
  p3_(idn,idn) = 1.;
  p3_(ivx,ivz) = 1.;
  p3_(ivy,ivx) = 1.;
  p3_(ivz,ivy) = 1.;
  p3_(ien,ien) = 1.;
  phydro->pmy_block->pdebug->Leave();
}

ImplicitSolver::~ImplicitSolver() {
  //delete[] usend_top_;
  //delete[] urecv_bot_;
  //delete[] usend_bot_;
  //delete[] urecv_top_;

  FreeCArray(buffer_);
  //FreeCArray(coefficients_);
  //FreeCArray(jacobian_);

#ifdef MPI_PARALLEL
  FreeCArray(req_send_data1_);
  FreeCArray(req_send_data2_);
  FreeCArray(req_send_data6_);
  FreeCArray(req_send_data7_);
#endif
}

void ImplicitSolver::SetDirection(CoordinateDirection dir) {
  mydir_ = dir;
  MeshBlock *pmb = pmy_hydro->pmy_block;
  int nc1, nc2, nc3;
  if (dir == X1DIR) {
    nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  } else if (dir == X2DIR) {
    nc1 = pmb->ncells2, nc2 = pmb->ncells3, nc3 = pmb->ncells1;
  } else { // X3DIR
    nc1 = pmb->ncells3, nc2 = pmb->ncells1, nc3 = pmb->ncells2;
  }

  du_.SetDim1(nc1);
  du_.SetDim2(nc2);
  du_.SetDim3(nc3);

  coefficients_.SetDim2(nc1);
  coefficients_.SetDim3(nc2);
  coefficients_.SetDim4(nc3);

  if ((pmb->pmy_mesh->mesh_bcs[2*dir] == BoundaryFlag::periodic) &&
     (pmb->pmy_mesh->mesh_bcs[2*dir+1] == BoundaryFlag::periodic)) {
    periodic_boundary = true;
  } else {
    periodic_boundary = false;
  }

  if (pmb->pbval->block_bcs[2*dir] == BoundaryFlag::polar)
    pole_at_bot = true;
  else
    pole_at_bot = false;

  if (pmb->pbval->block_bcs[2*dir+1] == BoundaryFlag::polar)
    pole_at_top = true;
  else
    pole_at_top = false;
}

void ImplicitSolver::FindNeighbors() {
  // find top and bot neighbor
  has_top_neighbor = false;
  has_bot_neighbor = false;
  first_block = true;
  last_block = true;

  for (int n = 0; n < pmy_hydro->pmy_block->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmy_hydro->pmy_block->pbval->neighbor[n];
    if (mydir_ == X1DIR) {
      if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
        bblock = nb;
        has_bot_neighbor = true;
      } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
        tblock = nb;
        has_top_neighbor = true;
      }
    } else if (mydir_ == X2DIR) {
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

  MeshBlock *pmb = pmy_hydro->pmy_block;
  if (mydir_ == X1DIR) {
    if (pmb->block_size.x1min > pmb->pmy_mesh->mesh_size.x1min)
      first_block = false;
    if (pmb->block_size.x1max < pmb->pmy_mesh->mesh_size.x1max)
      last_block = false;
  } else if (mydir_ == X2DIR) {
    if (pmb->block_size.x2min > pmb->pmy_mesh->mesh_size.x2min)
      first_block = false;
    if (pmb->block_size.x2max < pmb->pmy_mesh->mesh_size.x2max)
      last_block = false;
  } else { // X3DIR
    if (pmb->block_size.x3min > pmb->pmy_mesh->mesh_size.x3min)
      first_block = false;
    if (pmb->block_size.x3max < pmb->pmy_mesh->mesh_size.x3max)
      last_block = false;
  }

  //if (first_block)
  //  has_bot_neighbor = false;

  //if (last_block)
  //  has_top_neighbor = false;

  //if (pmb->pbval->block_bcs[2*mydir_] == BoundaryFlag::polar)
  //  first_block = true;
  //if (pmb->pbval->block_bcs[2*mydir_+1] == BoundaryFlag::polar)
  //  last_block = true;

  //std::cout << "dir = " << mydir_ << std::endl;
  //std::cout << "I'm rank " << Globals::my_rank << std::endl;
  //std::cout << "first_block " << first_block << std::endl;
  //std::cout << "last_block " << last_block << std::endl;
}

/*void ImplicitSolver::SynchronizeConserved(AthenaArray<Real> const& du,
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
      int tag = CreateMPITag(bblock.snb.gid, pmb->gid, "b");
      MPI_Isend(usend_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_sync_bot_);
      if (pmb->pbval->block_bcs[2*mydir_] == BoundaryFlag::polar)
        tag = CreateMPITag(pmb->gid, bblock.snb.gid, "b");
      else
        tag = CreateMPITag(pmb->gid, bblock.snb.gid, "t");
      MPI_Irecv(urecv_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_recv_sync_bot_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->phydro->pimp->urecv_top_, usend_bot_, sbot*sizeof(Real));
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
      int tag = CreateMPITag(tblock.snb.gid, pmb->gid, "t");
      MPI_Isend(usend_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_sync_top_);
      if (pmb->pbval->block_bcs[2*mydir_+1] == BoundaryFlag::polar)
        tag = CreateMPITag(pmb->gid, tblock.snb.gid, "t");
      else
        tag = CreateMPITag(pmb->gid, tblock.snb.gid, "b");
      MPI_Irecv(urecv_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_recv_sync_top_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->phydro->pimp->urecv_bot_, usend_top_, stop*sizeof(Real));
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
}*/

int ImplicitSolver::CreateMPITag(int recvid, int sendid, std::string phys) {
  //return (lid<<17) | (bufid<<11) | phys;
  std::string str = std::to_string(recvid);
  str += "x";
  str += std::to_string(sendid);
  str += "x";
  str += phys;
  return std::hash<std::string>{}(str)%(1<<24);
}

#undef MAX_DATA_SIZE
