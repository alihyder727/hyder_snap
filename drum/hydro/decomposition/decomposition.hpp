#ifndef DECOMPOSITION_HPP
#define DECOMPOSITION_HPP

// Athena++ headers
#include "../hydro.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

class Decomposition {
public:
// data
  Hydro *pmy_hydro;
  bool has_top_neighbor, has_bot_neighbor;
  NeighborBlock tblock, bblock;

// functions
  Decomposition(Hydro *phydro);
  ~Decomposition();
  void FindNeighbors();
  int CreateMPITag(int lid, int bufid, int phy);

  void RecvFromTop(AthenaArray<Real> &psf, int kl, int ku, int jl, int ju);
  void SendToBottom(AthenaArray<Real> const& psf, int kl, int ku, int jl, int ju);
  void SyncNewVariables(AthenaArray<Real> const& w, int kl, int ku, int jl, int ju);
  void WaitToFinishSync(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);

  void ChangeToBuoyancy(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);
  void RestoreFromBuoyancy(AthenaArray<Real> &w, AthenaArray<Real> &wl, AthenaArray<Real> &wr,
    int k, int j, int il, int iu);
  void ApplyBoundaryCondition(AthenaArray<Real> &w, AthenaArray<Real> &psf,
    int kl, int ku, int jl, int ju);

  void RecvFromTop(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
    int kl, int ku, int jl, int ju);
  void SendToBottom(AthenaArray<Real> &psf, AthenaArray<Real> &entropy,
    int kl, int ku, int jl, int ju);
  void WaitToFinishSend();

  void ChangeToPerturbation(AthenaArray<Real> &w, int kl, int ku, int jl, int ju);
  void RestoreFromPerturbation(AthenaArray<Real> &w, AthenaArray<Real> &wl, AthenaArray<Real> &wr,
    int k, int j, int il, int iu);

private:
  // pressure decomposition
  AthenaArray<Real> psf_;         // hydrostatic pressure at cell face
  AthenaArray<Real> pres_, dens_; // save of original w
  Real *buffer_;                  // MPI data buffer
  Real *wsend_top_, *wrecv_top_;  // MPI data buffer
  Real *wsend_bot_, *wrecv_bot_;  // MPI data buffer

  AthenaArray<Real> entropy_;     // pseudo entropy

#ifdef MPI_PARALLEL
  MPI_Request req_send_to_bot_;
  MPI_Request req_send_sync_top_;
  MPI_Request req_send_sync_bot_;
  MPI_Request req_recv_sync_top_;
  MPI_Request req_recv_sync_bot_;
#endif
};

#endif
