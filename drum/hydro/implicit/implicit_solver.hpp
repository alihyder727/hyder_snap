#ifndef IMPLICIT_SOLVER_HPP_
#define IMPLICIT_SOLVER_HPP_

// C/C++ headers
#include <vector>

// Athena++ headers
#include "../../athena.hpp"
#include "../../bvals/bvals_interfaces.hpp"

class ImplicitSolver {
  friend class Hydro;
public:
// data
  Hydro *pmy_hydro;
  bool has_top_neighbor, has_bot_neighbor;
  bool first_block, last_block;
  bool periodic_boundary;
  bool pole_at_bot, pole_at_top;
  NeighborBlock tblock, bblock;
  CoordinateDirection mydir;

// functions
  ImplicitSolver(Hydro *phydro, CoordinateDirection dir);
  ~ImplicitSolver();

// utility functions
  void FindNeighbors();
  int CreateMPITag(int lid, int bufid, std::string phy);

// correction methods
  void PartialCorrection(AthenaArray<Real>& du, AthenaArray<Real> const& w, Real dt);
  void FullCorrection(AthenaArray<Real>& du, AthenaArray<Real> const& w, Real dt);

// tri-diagonal solver
  template<typename T1, typename T2>
  void ForwardSweep(std::vector<T1> &a, std::vector<T1> &b, std::vector<T1> &c, 
    std::vector<T2> &delta, std::vector<T2> &corr, Real dt,
    int k, int j, int il, int iu);

  template<typename T1, typename T2>
  void BackwardSubstitution(std::vector<T1> &a, std::vector<T2> &delta, 
    int kl, int ku, int jl, int ju, int il, int iu);

// periodic solver
  template<typename T1, typename T2>
  void PeriodicForwardSweep(std::vector<T1> &a, std::vector<T1> &b, std::vector<T1> &c,
    std::vector<T2> &delta, std::vector<T2> &corr, Real dt,
    int k, int j, int il, int iu);

  template<typename T1, typename T2>
  void PeriodicBackwardSubstitution(std::vector<T1> &a, std::vector<T1> &c,
    std::vector<T2> &delta, int kl, int ku, int jl, int ju, int il, int iu);

// communications
  void SynchronizeConserved(AthenaArray<Real> const& du,
    int kl, int ku, int jl, int ju, int is, int ie);
  void WaitToFinishSync(int kl, int ku, int jl, int ju, int is, int ie);

  template<typename T>
  void SendBuffer(T const& a, int k, int j, NeighborBlock nb);

  template<typename T1, typename T2>
  void SendBuffer(T1 const& a, T2 const& b, int k, int j, NeighborBlock nb);

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void SendBuffer(T1 const& a, T2 const&b, T3 const& c, T4 const& d,
    T5 const& e, T6 const& f, int k, int j, NeighborBlock ntop);

  template<typename T1, typename T2, typename T3, 
           typename T4, typename T5, typename T6, typename T7>
  void SendBuffer(T1 const& a, T2 const&b, T3 const& c, T4 const& d,
    T5 const& e, T6 const& f, T7 const& g, int k, int j, NeighborBlock ntop);

  template<typename T>
  void RecvBuffer(T &a, int k, int j, NeighborBlock nb);

  template<typename T1, typename T2>
  void RecvBuffer(T1 &a, T2 &b, int k, int j, NeighborBlock nb);

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, 
    T5 &e, T6 &f, int k, int j, NeighborBlock nb);

  template<typename T1, typename T2, typename T3,
           typename T4, typename T5, typename T6, typename T7>
  void RecvBuffer(T1 &a, T2 &b, T3 &c, T4 &d, 
    T5 &e, T6 &f, T7 &g, int k, int j, NeighborBlock nb);

  template<typename T1, typename T2>
  void SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
    int k, int j, int il, int iu);

  template<typename T1, typename T2, typename T3, typename T4>
  void SaveCoefficients(std::vector<T1> &a, std::vector<T2> &b,
    std::vector<T3> &c, std::vector<T4> &d, int k, int j, int il, int iu);

  template<typename T1, typename T2>
  void LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
    int k, int j, int il, int iu);

  template<typename T1, typename T2, typename T3, typename T4>
  void LoadCoefficients(std::vector<T1> &a, std::vector<T2> &b,
    std::vector<T3> &c, std::vector<T4> &d, int k, int j, int il, int iu);

private:
  Real *usend_top_, *urecv_top_;  // MPI data buffer
  Real *usend_bot_, *urecv_bot_;  // MPI data buffer

  Real ***buffer_;                  // MPI data buffer
  Real ****coefficients_; // archive of coefficients in the tri-diagonal matrix
  AthenaArray<Real> du_;  // stores implicit solution

#ifdef MPI_PARALLEL
  MPI_Request **req_send_data1_;
  MPI_Request **req_send_data2_;
  MPI_Request **req_send_data6_;
  MPI_Request **req_send_data7_;
  MPI_Request req_send_sync_top_;
  MPI_Request req_send_sync_bot_;
  MPI_Request req_recv_sync_top_;
  MPI_Request req_recv_sync_bot_;
#endif
};

#endif
