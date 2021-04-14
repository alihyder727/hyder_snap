// C/C++ headers
#include <iostream>
#include <sstream>
#include <functional>

// Athena++ headers
#include "../../mesh/mesh.hpp"
#include "../../globals.hpp"
#include "decomposition.hpp"

Decomposition::Decomposition(Hydro *phydro):
  pmy_hydro(phydro), has_top_neighbor(false), has_bot_neighbor(false)
{
  MeshBlock *pmb = phydro->pmy_block;
  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;
  // allocate hydrostatic and nonhydrostatic pressure
  psf_.NewAthenaArray(nc3, nc2, nc1+1);
  pres_.NewAthenaArray(nc3, nc2, nc1);
  dens_.NewAthenaArray(nc3, nc2, nc1);

  buffer_ = new Real [(NGHOST+2)*nc3*nc2];
  wsend_top_ = new Real [2*NGHOST*nc3*nc2];
  wrecv_top_ = new Real [2*NGHOST*nc3*nc2];
  wsend_bot_ = new Real [2*NGHOST*nc3*nc2];
  wrecv_bot_ = new Real [2*NGHOST*nc3*nc2];

  // allocate polytropic index and pseudo entropy
  entropy_.NewAthenaArray(nc3, nc2);
}

Decomposition::~Decomposition()
{
  delete[] buffer_;
  delete[] wsend_top_;
  delete[] wrecv_top_;
  delete[] wsend_bot_;
  delete[] wrecv_bot_;
}

int Decomposition::CreateMPITag(int recvid, int sendid, int phys)
{
  std::string str = std::to_string(recvid);
  str += std::to_string(sendid);
  str += std::to_string(phys);
  return std::hash<std::string>{}(str)%(1<<24);
  //return (recvid<<11) | (sendid<<5) | phys;
}

void Decomposition::FindNeighbors()
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  has_top_neighbor = false;
  has_bot_neighbor = false;
  for (int n = 0; n < pmb->pbval->nneighbor; ++n) {
    NeighborBlock& nb = pmb->pbval->neighbor[n];
    if ((nb.ni.ox1 == -1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      bblock = nb;
      has_bot_neighbor = true;
    } if ((nb.ni.ox1 == 1) && (nb.ni.ox2 == 0) && (nb.ni.ox3 == 0)) {
      tblock = nb;
      has_top_neighbor = true;
    }
  }
}

// FIXME: local boundary has not been implemented
// Ordering the meshblocks need to be worked out such that
// the upper boundary executes before the lower boundary
void Decomposition::RecvFromTop(AthenaArray<Real> &psf, int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  int ssize = (ju-jl+1)*(ku-kl+1)*(NGHOST+1);

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (tblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmb->gid, tblock.snb.gid, 22);
    MPI_Recv(buffer_, ssize, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } else {  // local boundary
    // need to wait for the top boundary to finish
    msg << "### FATAL ERROR in Decompositin::RecvFromTop" << std::endl
        << "Local boundary not yet implemented" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = 0; i <= NGHOST; ++i)
        psf(k,j,pmb->ie+i+1) = buffer_[p++];
}

void Decomposition::SendToBottom(AthenaArray<Real> const& psf, int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  int ssize = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = 0; i <= NGHOST; ++i)
        buffer_[ssize++] = psf(k,j,pmb->is+i);

  if (bblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 22);
    MPI_Isend(buffer_, ssize, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_to_bot_);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->phydro->pdec->buffer_, buffer_, ssize*sizeof(Real));
  }
}

void Decomposition::SyncNewVariables(AthenaArray<Real> const& w, int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;

  if (has_bot_neighbor) {
    int sbot = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 0; i < NGHOST; ++i) {
          wsend_bot_[sbot++] = w(IDN,k,j,pmb->is+i);
          wsend_bot_[sbot++] = w(IPR,k,j,pmb->is+i);
        }
    if (bblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 17);
      MPI_Isend(wsend_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_sync_bot_);
      tag = CreateMPITag(pmb->gid, bblock.snb.gid, 19);
      MPI_Irecv(wrecv_bot_, sbot, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_recv_sync_bot_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->phydro->pdec->wrecv_top_, wsend_bot_, sbot*sizeof(Real));
      std::memcpy(wrecv_bot_, pbl->phydro->pdec->wsend_top_, sbot*sizeof(Real));
    }
  }

  if (has_top_neighbor) {
    int stop = 0;
    for (int k = kl; k <= ku; ++k)
      for (int j = jl; j <= ju; ++j)
        for (int i = 0; i < NGHOST; ++i) {
          wsend_top_[stop++] = w(IDN,k,j,pmb->ie-i);
          wsend_top_[stop++] = w(IPR,k,j,pmb->ie-i);
        }
    if (tblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(tblock.snb.gid, pmb->gid, 19);
      MPI_Isend(wsend_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_sync_top_);
      tag = CreateMPITag(pmb->gid, tblock.snb.gid, 17);
      MPI_Irecv(wrecv_top_, stop, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD,
        &req_recv_sync_top_);
#endif
    } else {  // local boundary
      MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
      std::memcpy(pbl->phydro->pdec->wrecv_bot_, wsend_top_, stop*sizeof(Real));
      std::memcpy(wrecv_top_, pbl->phydro->pdec->wsend_bot_, stop*sizeof(Real));
    }
  }
}

void Decomposition::WaitToFinishSend()
{
#ifdef MPI_PARALLEL
  MPI_Status status;
  if (has_bot_neighbor && (bblock.snb.rank != Globals::my_rank))
    MPI_Wait(&req_send_to_bot_, &status);
#endif
}

void Decomposition::WaitToFinishSync(AthenaArray<Real> &w,
  int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
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
        for (int i = 1; i <= NGHOST; ++i) {
          w(IDN,k,j,pmb->is-i) = wrecv_bot_[p++];
          w(IPR,k,j,pmb->is-i) = wrecv_bot_[p++];
        }
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
        for (int i = 1; i <= NGHOST; ++i) {
          w(IDN,k,j,pmb->ie+i) = wrecv_top_[p++];
          w(IPR,k,j,pmb->ie+i) = wrecv_top_[p++];
        }
  }
}

// FIXME: local boundary has not been implemented
// Ordering the meshblocks need to be worked out such that
// the upper boundary executes before the lower boundary
void Decomposition::RecvFromTop(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
  int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  int ssize = (NGHOST+2)*(ju-jl+1)*(ku-kl+1);

  std::stringstream msg;
#ifdef MPI_PARALLEL
  MPI_Status status;
#endif

  if (tblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(pmb->gid, tblock.snb.gid, 22);
    MPI_Recv(buffer_, ssize, MPI_ATHENA_REAL, tblock.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
  } else {  // local boundary
    // need to wait for the top boundary to finish
    msg << "### FATAL ERROR in Decompositin::RecvFromTop" << std::endl
        << "Local boundary not yet implemented" << std::endl;
    ATHENA_ERROR(msg);
  }
  int p = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      for (int i = 0; i <= NGHOST; ++i)
        psf(k,j,pmb->ie+i+1) = buffer_[p++];
      //psf(k,j,pmb->ie+1) = buffer_[p++];
      entropy(k,j) = buffer_[p++];
    }
}

void Decomposition::SendToBottom(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
  int kl, int ku, int jl, int ju)
{
  MeshBlock *pmb = pmy_hydro->pmy_block;
  int ssize = 0;
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j) {
      for (int i = 0; i <= NGHOST; ++i)
        buffer_[ssize++] = psf(k,j,pmb->is+i);
      //buffer_[ssize++] = psf(k,j,pmb->is);
      buffer_[ssize++] = entropy(k,j);
    }

  if (bblock.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
    int tag = CreateMPITag(bblock.snb.gid, pmb->gid, 22);
    MPI_Isend(buffer_, ssize, MPI_ATHENA_REAL, bblock.snb.rank, tag, MPI_COMM_WORLD,
      &req_send_to_bot_);
#endif
  } else {  // local boundary
    MeshBlock *pbl = pmy_hydro->pmy_block->pmy_mesh->FindMeshBlock(bblock.snb.gid);
    std::memcpy(pbl->phydro->pdec->buffer_, buffer_, ssize*sizeof(Real));
  }
}
