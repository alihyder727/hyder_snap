#ifndef IMPLICIT_SOLVER_HPP_
#define IMPLICIT_SOLVER_HPP_

// C/C++ headers
#include <cstring>

// Athena++ headers
#include "../../athena.hpp"
#include "../../globals.hpp"
#include "../../bvals/bvals_interfaces.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"

// MPI headers
#ifdef MPI_PARALLEL
  #include <mpi.h>
#endif

class ImplicitSolver {
  friend class Hydro;
public:
// data
  Hydro *pmy_hydro;
  bool has_top_neighbor, has_bot_neighbor;
  NeighborBlock tblock, bblock;
  CoordinateDirection mydir;

// functions
  ImplicitSolver(Hydro *phydro, CoordinateDirection dir);
  ~ImplicitSolver();
  void FindNeighbors();
  int CreateMPITag(int lid, int bufid, int phy);

  void WaitSendTop(int kl, int ku, int jl, int ju);
  void WaitSendBot(int kl, int ku, int jl, int ju);
  void WaitToFinishSend(int kl, int ku, int jl, int ju);

  void SynchronizeConserved(AthenaArray<Real> const& du,
    int kl, int ku, int jl, int ju, int is, int ie);
  void WaitToFinishSync(int kl, int ku, int jl, int ju, int is, int ie);

  void PartialCorrection(AthenaArray<Real> const& du, AthenaArray<Real> const& w, Real dt);

  template<typename T1, typename T2>
  void ForwardSweep(std::vector<T1> &a, std::vector<T1> &b, std::vector<T1> &c, 
    std::vector<T2> &delta, std::vector<T2> &corr, Real dt,
    int k, int j, int il, int iu);

  template<typename T1, typename T2>
  void BackwardSubstitution(std::vector<T1> &a, std::vector<T2> &delta, 
    int kl, int ku, int jl, int ju, int il, int iu);

  template<typename T>
  void SendBotBuffer(T &a, int k, int j, NeighborBlock nbot) {
    int s1 = a.size();
    int phy = k << 7 | j << 1 | 0;

    memcpy(buffer_, a.data(), s1*sizeof(Real));

    if (nbot.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(nbot.snb.gid, pmy_hydro->pmy_block->gid, phy);
      MPI_Isend(buffer_, s1, MPI_ATHENA_REAL, nbot.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_bot_data_[k][j]);
#endif
    } // local boundary
  }

  template<typename T>
  void RecvTopBuffer(T &a, int k, int j, NeighborBlock ntop) {
    int s1 = a.size();
    int phy = k << 7 | j << 1 | 0;
#ifdef MPI_PARALLEL
    MPI_Status status;
#endif

    if (ntop.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(pmy_hydro->pmy_block->gid, ntop.snb.gid, phy);
      MPI_Recv(buffer_, s1, MPI_ATHENA_REAL, ntop.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
    } // local boundary

    memcpy(a.data(), buffer_, s1*sizeof(Real));
  }

  template<typename T1, typename T2>
  void SendTopBuffer(T1 &a, T2 &b, int k, int j, NeighborBlock ntop) {
    int s1 = a.size(), s2 = b.size();
    int phy = k << 7 | j << 1 | 1;

    memcpy(buffer_, a.data(), s1*sizeof(Real));
    memcpy(buffer_ + s1, b.data(), s2*sizeof(Real));

    if (ntop.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(ntop.snb.gid, pmy_hydro->pmy_block->gid, phy);
      MPI_Isend(buffer_, s1+s2, MPI_ATHENA_REAL, ntop.snb.rank, tag, MPI_COMM_WORLD,
        &req_send_top_data_[k][j]);
#endif
    } // local boundary
  }

  template<typename T1, typename T2>
  void RecvBotBuffer(T1 &a, T2 &b, int k, int j, NeighborBlock nbot) {
    int s1 = a.size(), s2 = b.size();
    int phy = k << 7 | j << 1 | 1;
#ifdef MPI_PARALLEL
    MPI_Status status;
#endif

    if (nbot.snb.rank != Globals::my_rank) { // MPI boundary
#ifdef MPI_PARALLEL
      int tag = CreateMPITag(pmy_hydro->pmy_block->gid, nbot.snb.gid, phy);
      MPI_Recv(buffer_, s1+s2, MPI_ATHENA_REAL, nbot.snb.rank, tag, MPI_COMM_WORLD, &status);
#endif
    } // local boundary

    memcpy(a.data(), buffer_, s1*sizeof(Real));
    memcpy(b.data(), buffer_ + s1, s2*sizeof(Real));
  }

  template<typename T1, typename T2>
  void SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
    int k, int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      int s1 = a[i].size(), s2 = b[i].size();
      memcpy(coefficients_[k][j][i], a[i].data(), s1*sizeof(Real));
      memcpy(coefficients_[k][j][i] + s1, b[i].data(), s2*sizeof(Real));
    }
  }

  template<typename T1, typename T2>
  void LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
    int k, int j, int il, int iu) {
    for (int i = il; i <= iu; ++i) {
      int s1 = a[i].size(), s2 = b[i].size();
      memcpy(a[i].data(), coefficients_[k][j][i], s1*sizeof(Real));
      memcpy(b[i].data(), coefficients_[k][j][i] + s1, s2*sizeof(Real));
    }
  }

private:
  Real *buffer_;                  // MPI data buffer
  Real *usend_top_, *urecv_top_;  // MPI data buffer
  Real *usend_bot_, *urecv_bot_;  // MPI data buffer
  Real ****coefficients_; // archive of coefficients in the tri-diagonal matrix
  AthenaArray<Real> du_;  // stores implicit solution

#ifdef MPI_PARALLEL
  MPI_Request **req_send_bot_data_;
  MPI_Request **req_send_top_data_;
  MPI_Request req_send_sync_top_;
  MPI_Request req_send_sync_bot_;
  MPI_Request req_recv_sync_top_;
  MPI_Request req_recv_sync_bot_;
#endif
};

#endif
